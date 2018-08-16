#ifndef SPACEX_H
#define SPACEX_H
#include "solver.h"
#include "particles.h"
#include "kahanNumber.h"
#include "macros.h"


#include "tools/timer.h"

#include "particleSystem.h"
#include "particleSystem/reguSystem/regularization.h"
#include "particleSystem/reguSystem/reguSystem.h"
#include "particleSystem/reguSystem/regularParticles.h"

#include "particleSystem/ARChain/ARchainSystem.h"
#include "particleSystem/ARChain/chain.h"
#include "particleSystem/ARChain/chainParticles.h"

#include "integrator/symplectic/symplectic2th.h"
#include "integrator/symplectic/symplectic4th.h"
#include "integrator/symplectic/symplectic6th.h"
#include "integrator/symplectic/symplectic8th.h"
#include "integrator/symplectic/symplectic10th.h"

#include "integrator/Gauss-Dadau.h"

#include "ODEiterator/BSIterator_impl3.h"
#include "ODEiterator/constIterator.h"
#include "ODEiterator/IAS15.h"

#include "interaction/interaction.h"
#include "interaction/force.h"
#include "interaction/postNewtonian.h"

//#include "wrapper/singleRun.h"
namespace SpaceH
{
    
    template
    <
        typename BasicF,
        typename VelForce = void,
        typename ExtPosForce = void,
        typename ExtVelForce = void,
        template<typename> class Regularitor = SpaceH::LogH
    >
    using ARchain = SpaceH::ARchainSystem
    <
        SpaceH::ChainParticles
        <
            typename BasicF::type::Scalar,
            BasicF::type::arraySize,
            !std::is_void<VelForce>::value || !std::is_void<ExtVelForce>::value
        >,
    
        SpaceH::Interactions
        <
            BasicF, VelForce, ExtPosForce, ExtVelForce
        >,
    
        Regularitor
        <
            SpaceH::ChainParticles
            <
                typename BasicF::type::Scalar,
                BasicF::type::arraySize,
                !std::is_void<VelForce>::value || !std::is_void<ExtVelForce>::value
            >
        >
    >;
    
    template
    <
        typename BasicF,
        typename VelForce = void,
        typename ExtPosForce = void,
        typename ExtVelForce = void,
        template<typename> class Regularitor = SpaceH::LogH
    >
    using GAR = SpaceH::ReguSystem
    <
        SpaceH::ReguParticles
        <
            typename BasicF::type::Scalar,
            BasicF::type::arraySize,
            !std::is_void<VelForce>::value || !std::is_void<ExtVelForce>::value
        >,
    
        SpaceH::Interactions
        <
            BasicF, VelForce, ExtPosForce, ExtVelForce
        >,
    
        Regularitor
        <
            SpaceH::ReguParticles
            <
                typename BasicF::type::Scalar,
                BasicF::type::arraySize,
                !std::is_void<VelForce>::value || !std::is_void<ExtVelForce>::value
            >
        >
    >;
    
    template
    <
        typename BasicF,
        typename VelForce = void,
        typename ExtPosForce = void,
        typename ExtVelForce = void
    >
    using Basic = SpaceH::ParticleSystem
    <
        SpaceH::Particles
        <
            typename BasicF::type::Scalar,
            BasicF::type::arraySize,
            !std::is_void<VelForce>::value || !std::is_void<ExtVelForce>::value
        >,
    
        SpaceH::Interactions
        <
            BasicF, VelForce, ExtPosForce, ExtVelForce
        >
    >;
    
    /**  @brief Alias of template name, linking the particle system, integrator and ODE iterator*/
    template
    <
        typename                           ParticSys,
        template<typename, typename> class ODEiterator = SpaceH::BSIterator,
        template<typename> class           Integrator  = SpaceH::symplectic2th
    >
    using Nbody = SpaceH::Solver
    <
        ParticSys,
    
        ODEiterator
        <
            ParticSys,
            Integrator<ParticSys>
        >
    >;
}
#endif
