
#ifndef GAUSSRADAU_H
#define GAUSSRADAU_H
#include "PredictorCorrector.h"

namespace SpaceH
{
    namespace RADAU
    {
        /** Gauss radau stepping. see details in
         https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
         */
        constexpr double H[8] = {0, 0.05626256053692215, 0.18024069173689236, 0.35262471711316964, 0.54715362633055538, 0.73421017721541053, 0.88532094683909577, 0.97752061356128750};
        
        constexpr double c[21] =
        {
            -0.05626256053692215,  0.01014080283006363, -0.23650325227381452, -0.00357589772925162,  0.09353769525946207,
            -0.58912796938698420,  0.00195656540994722, -0.05475538688906869,  0.41588120008230689, -1.13628159571753962,
            -0.00143653023637089,  0.04215852772126871, -0.36009959650205681,  1.25015071184069093, -1.87049177293294999,
             0.00127179030902687, -0.03876035791590677,  0.36096224345284600, -1.46688420840042699,  2.90613625930842900,
            -2.75581271977204567
        };
        
        constexpr double d[21] =
        {
            0.05626256053692215, 0.00316547571817083, 0.23650325227381452, 0.00017809776922174, 0.04579298550602792,
            0.58912796938698420, 0.00001002023652233, 0.00843185715352570, 0.25353406905456927, 1.13628159571753962,
            0.00000056376416393, 0.00152978400250047, 0.09783423653244401, 0.87525466468409108, 1.87049177293294999,
            0.00000003171881540, 0.00027629309098265, 0.03602855398373646, 0.57673300027707874, 2.24858876076915948,
            2.75581271977204567
        };
        
        constexpr double r[28] =
        {
           17.77380891407799979, 5.54813671853721679, 8.06593864838188779, 2.83587607864443880, 3.37424997696263551,
            5.80100155926406202, 1.82764026751759778, 2.03711183535858442, 2.72544221180822621, 5.14062410581093232,
            1.36200781606246957, 1.47504021756041159, 1.80515358014025140, 2.62064492638703506, 5.34597689987110947,
            1.12953387533678984, 1.20618766605844563, 1.41827826373473909, 1.87724249618681016, 2.95711601729045581,
            6.61766201370242158, 1.02299632982348676, 1.08547219393864247, 1.25426462228187785, 1.60026654949081615,
            2.32359830021969449, 4.10997577834455807, 10.84602619023684689
        };
        
        
    }
    
    /** @brief Gauss RADAU iterator
     *
     */
    template <typename ParticSys, typename Integrator>
    class GaussRADAU : public PCIterator
    {
    public:
        /* Typedef */
        using type = typename ParticSys::type
        using Scalar = typename type::Scalar;
        using PlainArray =  typename type::PlainArray;
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
        Scalar iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength);
        
    private:
        
        ScalarBuffer
    };

}
#endif
