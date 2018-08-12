
#ifndef GAUSSDADAU_H
#define GAUSSDADAU_H
#include "../devTools.h"
namespace SpaceH
{
    namespace Radau
    {
        /** Gauss radau stepping. see details in
         https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
         */
        
        constexpr double interplt[8] = {5.62625605369221465e-02, 1.80240691736892365e-01, 3.52624717113169637e-01,   5.47153626330555383e-01, 7.34210177215410532e-01, 8.85320946839095768e-01, 9.77520613561287502e-01};
        
        constexpr double VC[8][8] =
        {
            {5.62625605369221465e-02, 1.58273785908541486e-03, 5.93659230739144767e-05, 2.50505913058228235e-06, 1.12752832786364176e-07, 5.28646923362689616e-09, 2.54940253100151162e-10, 1.25506424954273203e-11},
            {1.80240691736892365e-01, 1.62433534788967299e-02, 1.95180884477546917e-03, 2.63846532240386510e-04, 3.80447051867100280e-05, 5.71433664981562659e-06, 8.82819420497352521e-07, 1.39229985150554609e-07},
            {3.52624717113169637e-01, 6.21720955595714483e-02, 1.46156117393512197e-02, 3.86536946626848347e-03, 1.09041985166464595e-03, 3.20424159773191781e-04, 9.68481245967829830e-05, 2.98821622215214056e-05},
            {5.47153626330555383e-01, 1.49688545403338535e-01, 5.46017536250551147e-02, 2.24066606249673402e-02, 9.80790849192715933e-03, 4.47202724839682909e-03, 2.09733079372232626e-03, 1.00411688072492376e-03},
            {7.34210177215410532e-01, 2.69532292163342235e-01, 1.31928901329682186e-01, 7.26476565188252793e-02, 4.26709190135767871e-02, 2.61078525090855350e-02, 1.64302723006367057e-02, 1.05553639953544359e-02},
            {8.85320946839095768e-01, 3.91896589456036559e-01, 2.31302839760153783e-01, 1.53582936827273253e-01, 1.08776152840200480e-01, 8.02515055216670714e-02, 6.08985761603187542e-02, 4.71754369689804109e-02},
            {9.77520613561287502e-01, 4.77773274968617989e-01, 3.11355483260339461e-01, 2.28267302274238665e-01, 1.78508794700074913e-01, 1.45413355434419272e-01, 1.21838187792222094e-01, 1.04211922575117272e-01},
            {1,                       5.00000000000000000e-01, 3.33333333333333315e-01, 2.50000000000000000e-01, 2.00000000000000011e-01, 1.66666666666666657e-01, 1.42857142857142849e-01, 1.25000000000000000e-01}
        };
        
        constexpr double PC[8][9] =
        {
            {5.62625605369221465e-02, 1.58273785908541486e-03, 2.96829615369572383e-05, 8.35019710194094082e-07, 2.81882081965910440e-08, 1.05729384672537927e-09, 4.24900421833585270e-11, 1.79294892791818853e-12, 7.84590314640274775e-14},
            {1.80240691736892365e-01, 1.62433534788967299e-02, 9.75904422387734584e-04, 8.79488440801288367e-05, 9.51117629667750701e-06, 1.14286732996312544e-06, 1.47136570082892105e-07, 1.98899978786506594e-08, 2.78832320378369031e-09},
            {3.52624717113169637e-01, 6.21720955595714483e-02, 7.30780586967560986e-03, 1.28845648875616108e-03, 2.72604962916161487e-04, 6.40848319546383509e-05, 1.61413540994638316e-05, 4.26888031736020116e-06, 1.17079877778820332e-06},
            {5.47153626330555383e-01, 1.49688545403338535e-01, 2.73008768125275574e-02, 7.46888687498911338e-03, 2.45197712298178983e-03, 8.94405449679365819e-04, 3.49555132287054341e-04, 1.43445268674989096e-04, 6.10451325053741953e-05},
            {7.34210177215410532e-01, 2.69532292163342235e-01, 6.59644506648410928e-02, 2.42158855062750943e-02, 1.06677297533941968e-02, 5.22157050181710700e-03, 2.73837871677278443e-03, 1.50790914219349091e-03, 8.61095074400260456e-04},
            {8.85320946839095768e-01, 3.91896589456036559e-01, 1.15651419880076892e-01, 5.11943122757577487e-02, 2.71940382100501199e-02, 1.60503011043334129e-02, 1.01497626933864590e-02, 6.73934813842577262e-03, 4.64060028054731274e-03},
            {9.77520613561287502e-01, 4.77773274968617989e-01, 1.55677741630169730e-01, 7.60891007580795503e-02, 4.46271986750187283e-02, 2.90826710868838552e-02, 2.03063646320370168e-02, 1.48874175107310391e-02, 1.13188113884477790e-02},
            {1,                       5.00000000000000000e-01, 1.66666666666666657e-01, 8.33333333333333287e-02, 5.00000000000000028e-02, 3.33333333333333329e-02, 2.38095238095238082e-02, 1.78571428571428562e-02, 1.38888888888888881e-02}
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
        
        template<typename RadauTab>
        void transB2G(RadauTab& B, RadauTab& G)
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
        void transG2B(RadauTab& G, RadauTab& B)
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
    }
    
    /** @brief Gauss-Dadau integrator */
    template <typename ParticSys>
    class GaussDadau
    {
    public:
        /* Typedef */
        using type        = typename ParticSys::type;
        using Scalar      = typename type::Scalar;
        using Vector      = typename type::Vector;
        using VectorArray = typename type::VectorArray;
        
        template<typename T, size_t S>
        using Container = typename type::template Container<T, S>;
        
        using RadauArray = Container<Vector,7>;
        using RadauTab   = Container<RadauArray, ParticSys::arraySize>;
        /* Typedef */
        
        /*Template parameter check*/
        /*Template parameter check*/
        
        /** @brief Order of the integrator*/
        static const int order{15};
        static const size_t finalPoint{7};
        
        /** @brief Interface to integrate particle system
         *
         *  This function integrate the particle system for one step with Gauss-Radau stepping.
         *  @param particles  Particle system need to be integrated.
         *  @param stepLength Step size for integration.
         */
        void integrate(ParticSys& particles, Scalar stepLength)
        {
            checkTabVolume<ParticSys::arraySize>(particles.particleNumber());
            calcuBTab(particles, stepLength);
            evaluateSystemAt(particles, stepLength, finalPoint);
        }
        
        void calcuBTab(const ParticSys& particles, Scalar stepLength)
        {
            for(size_t i = 0 ; i < finalPoint; ++i)
            {
                localSystem_ = particles;
                evaluateSystemAt(localSystem_, stepLength, i);
                calcuGTab(particles.acc(), localSystem_.acc(), i);
            }
            Radau::transG2B(Gtab_, Btab_);
        }
        
        void evaluateSystemAt(ParticSys& particleSys, Scalar stepLength, size_t index)
        {
            VectorArray dpos;
            VectorArray dvel;
            evaluateVelIncrement(dvel, particleSys.acc(), index);
            evaluatePosIncrement(dpos, particleSys.vel(), particleSys.acc(), stepLength, index);
            particleSys.advanceVel(dvel, stepLength);
            particleSys.advancePos(dpos, stepLength);
            particleSys.evaluateAcc();
        }
        
        inline void setBTab(RadauTab& iterBTab)
        {
            Btab_ = iterBTab;
        }
        
        inline const RadauTab& getBTab() const
        {
            return Btab_;
        }
        
        template<size_t isDYNAMICAL>
        typename std::enable_if<isDYNAMICAL != SpaceH::DYNAMICAL>::type
        checkTabVolume(size_t particleNum){}
        
        template<size_t isDYNAMICAL>
        typename std::enable_if<isDYNAMICAL == SpaceH::DYNAMICAL>::type
        checkTabVolume(size_t particleNum)
        {
            if(particleNum > particleNumber_)
            {
                Btab_.resize(particleNum);
                Gtab_.resize(particleNum);
                for(size_t i = 0 ; i < particleNum; ++i)//once particle number changes, old b Value should be droped away.
                {
                    Btab_[i][0].setZero();
                    Btab_[i][1].setZero();
                    Btab_[i][2].setZero();
                    Btab_[i][3].setZero();
                    Btab_[i][4].setZero();
                    Btab_[i][5].setZero();
                    Btab_[i][6].setZero();
                }
                particleNumber_ = particleNum;
            }
        }
    private:
        void evaluatePosIncrement(VectorArray& dpos, const VectorArray& vel, const VectorArray& acc, Scalar stepLength, size_t iter)
        {
            using namespace Radau;
            size_t size = vel.size();
            dpos.resize(size);
            
            for(size_t i = 0 ; i < size; ++i)
            {
                dpos[i] = vel[i]   *  PC[iter][0] + (acc[i]  *  PC[iter][1] + Btab_[i][0]*PC[iter][2] + Btab_[i][1]*PC[iter][3] + Btab_[i][2]*PC[iter][4]
                        + Btab_[i][3]*PC[iter][5] + Btab_[i][4]*PC[iter][6] + Btab_[i][5]*PC[iter][7] + Btab_[i][6]*PC[iter][8])*stepLength;
            }
        }
        
        void evaluateVelIncrement(VectorArray& dvel, const VectorArray& acc, size_t iter)
        {
            using namespace Radau;
            size_t size = acc.size();
            dvel.resize(size);
            for(size_t i = 0 ; i < size; ++i)
            {
                dvel[i] = acc[i]   *  VC[iter][0] + Btab_[i][0]*VC[iter][1] + Btab_[i][1]*VC[iter][2] + Btab_[i][2]*VC[iter][3]
                        + Btab_[i][3]*VC[iter][4] + Btab_[i][4]*VC[iter][5] + Btab_[i][5]*VC[iter][6] + Btab_[i][6]*VC[iter][7];
            }
        }
        
        void calcuGTab(const VectorArray& a0, const VectorArray& a, size_t iter)
        {
            size_t size   = a0.size();
            size_t offset = iter*(iter+1)/2 + 1;
            for(size_t i = 0 ; i < size; ++i)
            {
                Gtab_[i][iter] = (a[i] - a0[i])*Radau::gg[offset - 1];
                
                for(size_t j = 0 ; j < iter; ++j)
                    Gtab_[i][iter] -= Gtab_[i][j]*Radau::gg[j + offset];
            }
        }
        
    private:
        RadauTab Btab_;
        RadauTab Gtab_;
        ParticSys localSystem_;
        size_t particleNumber_{ParticSys::arraySize};
    };
    
}
#endif
