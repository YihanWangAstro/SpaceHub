
#ifndef CONSTITERATOR_H
#define CONSTITERATOR_H
namespace SpaceH
{
    /** @brief Most common iterator
     *
     *  Constant iterator keep the step length constant and integrate the particle system for one step.
     */
    template <typename ParticSys, typename Integrator>
    class constIterator
    {
    public:
        /* Typedef */
        using type   = typename ParticSys::type;
        using Scalar = typename type::Scalar;
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
            integrator.integrate(particles, stepLength);
            return stepLength;
        }
    private:
        Integrator integrator;
    };
}
#endif
