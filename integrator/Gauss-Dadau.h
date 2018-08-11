
#ifndef IAS15_H
#define IAS15_H
#include "../../devTools.h"
namespace SpaceH
{
    
    constexpr double VC[7][7] =
    {
        {1.58273785908541486e-03, 5.93659230739144767e-05, 2.50505913058228235e-06, 1.12752832786364176e-07, 5.28646923362689616e-09, 2.54940253100151162e-10, 1.25506424954273203e-11},
        {1.62433534788967299e-02, 1.95180884477546917e-03, 2.63846532240386510e-04, 3.80447051867100280e-05, 5.71433664981562659e-06, 8.82819420497352521e-07, 1.39229985150554609e-07},
        {6.21720955595714483e-02, 1.46156117393512197e-02, 3.86536946626848347e-03, 1.09041985166464595e-03, 3.20424159773191781e-04, 9.68481245967829830e-05, 2.98821622215214056e-05},
        {1.49688545403338535e-01, 5.46017536250551147e-02, 2.24066606249673402e-02, 9.80790849192715933e-03, 4.47202724839682909e-03, 2.09733079372232626e-03, 1.00411688072492376e-03},
        {2.69532292163342235e-01, 1.31928901329682186e-01, 7.26476565188252793e-02, 4.26709190135767871e-02, 2.61078525090855350e-02, 1.64302723006367057e-02, 1.05553639953544359e-02},
        {3.91896589456036559e-01, 2.31302839760153783e-01, 1.53582936827273253e-01, 1.08776152840200480e-01, 8.02515055216670714e-02, 6.08985761603187542e-02, 4.71754369689804109e-02},
        {4.77773274968617989e-01, 3.11355483260339461e-01, 2.28267302274238665e-01, 1.78508794700074913e-01, 1.45413355434419272e-01, 1.21838187792222094e-01, 1.04211922575117272e-01}
    }
    
    constexpr double PC[7][7] =
    {
        {2.96829615369572383e-05, 8.35019710194094082e-07, 2.81882081965910440e-08, 1.05729384672537927e-09, 4.24900421833585270e-11, 1.79294892791818853e-12, 7.84590314640274775e-14},
        {9.75904422387734584e-04, 8.79488440801288367e-05, 9.51117629667750701e-06, 1.14286732996312544e-06, 1.47136570082892105e-07, 1.98899978786506594e-08, 2.78832320378369031e-09},
        {7.30780586967560986e-03, 1.28845648875616108e-03, 2.72604962916161487e-04, 6.40848319546383509e-05, 1.61413540994638316e-05, 4.26888031736020116e-06, 1.17079877778820332e-06},
        {2.73008768125275574e-02, 7.46888687498911338e-03, 2.45197712298178983e-03, 8.94405449679365819e-04, 3.49555132287054341e-04, 1.43445268674989096e-04, 6.10451325053741953e-05},
        {6.59644506648410928e-02, 2.42158855062750943e-02, 1.06677297533941968e-02, 5.22157050181710700e-03, 2.73837871677278443e-03, 1.50790914219349091e-03, 8.61095074400260456e-04},
        {1.15651419880076892e-01, 5.11943122757577487e-02, 2.71940382100501199e-02, 1.60503011043334129e-02, 1.01497626933864590e-02, 6.73934813842577262e-03, 4.64060028054731274e-03},
        {1.55677741630169730e-01, 7.60891007580795503e-02, 4.46271986750187283e-02, 2.90826710868838552e-02, 2.03063646320370168e-02, 1.48874175107310391e-02, 1.13188113884477790e-02}
    }
    
    double velTaylorCoef(double h, size_t index)//make it const later
    {
        double var = 1;
        size_t order = index + 2;
        for(size_t i = 0 ; i < order; ++i)
            var*=h;
        
        return var/order;
    }
    
    double posTaylorCoef(double h, size_t index)//make it const later
    {
        double var = 1;
        size_t order = index + 3
        for(size_t i = 0 ; i < order; ++i)
            var*=h;
        
        return var/(order*(order-1));
    }
    
    void calcuVel(VectorArray& vel0, VectorArray& acc0 Scalar StepSize, size_t step)
    {
        size_t size = vel0.size();
        Scalar h = Radau::H[step];
        for(size_t i = 0 ; i < size; ++i)
        {
            this->eva_vel[i] = vel0[i] + (acc0[i] * h + cofBtab_[i][0]*VC[i][0] + cofBtab_[i][0]*VC[i][0] + cofBtab_[i][0]*VC[i][0]
                                          + cofBtab_[i][0]*VC[i][0] + cofBtab_[i][0]*VC[i][0] + cofBtab_[i][0]*VC[i][0] + cofBtab_[i][0]*VC[i][0])*stepSize;
        }
    }
    
    void calcuPos(VectorArray& pos0, VectorArray& vel0, VectorArray& acc0, Scalar StepSize, size_t step)
    {
        size_t size = pos0.size();
        Scalar h = Radau::H[step];
        for(size_t i = 0 ; i < size; ++i)
        {
            this->eva_pos[i] = pos0[i] + vel0[i]*h*stepSize + (acc0[i]*h*h*0.5 + cofBtab_[i][0]*PC[i][0] + cofBtab_[i][0]*PC[i][0] + cofBtab_[i][0]*PC[i][0]
                                                               + cofBtab_[i][0]*PC[i][0] + cofBtab_[i][0]*PC[i][0] + cofBtab_[i][0]*PC[i][0] + cofBtab_[i][0]*PC[i][0])*stepSize*stepSize;
        }
    }
    
    template<typename Dtype, size_t ArraySize>
    class DadauTab
    {
    public:
        /* Typedef */
        using type = SpaceH::ProtoType<Dtype, ArraySize>;
        using Vector = typename type::Vector;
        
        template<typename T, size_t S>
        using Container = typename type::template Container<T, S>;
        
        using RadauArray = Container<Vector,7>;
        /* Typedef */
        Contanier<RadauArray, ArraySize>;
    };
    /** @brief Second order symplectic integrator */
    template <typename ParticSys>
    class Dadau
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
        using RadauTab = Contanier<RadauArray, ParticSys::arraySize>;
        /* Typedef */
        
        /*Template parameter check*/
        /*Template parameter check*/
        
        /** @brief Order of the integrator*/
        static const int order{15};
        void integrate(ParticSys& particles, Scalar stepLength);
    private:
        RadauTab* Btab_;
    };
    
    /** @brief Interface to integrate particle system
     *
     *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
     *  @param particles  Particle system need to be integrated.
     *  @param stepLength Step size for integration.
     */
    template <typename ParticSys>
    void Dadau<ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
    {
        particles.drift(0.5 * stepLength);
        particles.kick(stepLength);
        particles.drift(0.5 * stepLength);
    }

}
#endif
