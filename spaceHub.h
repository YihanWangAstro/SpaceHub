#ifndef SPACEX_H
#define SPACEX_H

#include "solver.h"
#include "particles.h"
#include "kahanNumber.h"
#include "macros.h"

#include "tools/timer.h"
#include "callbacks.h"

#include "particleSystem.h"
#include "particleSystem/regularization.h"
#include "particleSystem/reguSystem.h"
#include "particleSystem/GARSystem.h"

#include "particleSystem/chain.h"
#include "particleSystem/chainParticles.h"

#include "integrator/symplectic/symplectic2th.h"
#include "integrator/symplectic/symplectic4th.h"
#include "integrator/symplectic/symplectic6th.h"
#include "integrator/symplectic/symplectic8th.h"
#include "integrator/symplectic/symplectic10th.h"

#include "integrator/Gauss-Dadau.h"

#include "ODEiterator/BSIterator_impl0.h"
#include "ODEiterator/constIterator.h"
#include "ODEiterator/IAS15.h"

#include "interaction/interaction.h"
#include "interaction/force.h"
#include "interaction/postNewtonian.h"

//#include "wrapper/singleRun.h"
namespace SpaceH {

    using DEFAULT_TYPE_CLASS = SpaceH::ProtoType<kahan<double>, SpaceH::DYNAMICAL>;

    template<
            typename BasicF      = SpaceH::NewtonianForce<DEFAULT_TYPE_CLASS>,
            typename VelForce    = void,
            typename ExtPosForce = void,
            typename ExtVelForce = void,
            template<typename> class Regularitor = SpaceH::LogH
    >
    using ARchain = SpaceH::GARSystem<
            SpaceH::ChainParticles<typename BasicF::type>,

            SpaceH::Interactions<BasicF, VelForce, ExtPosForce, ExtVelForce>,

            Regularitor<typename BasicF::type>
    >;

    template<
            typename BasicF      = SpaceH::NewtonianForce<DEFAULT_TYPE_CLASS>,
            typename VelForce    = void,
            typename ExtPosForce = void,
            typename ExtVelForce = void,
            template<typename> class Regularitor = SpaceH::LogH
    >
    using GAR = SpaceH::GARSystem<
            SpaceH::Particles<typename BasicF::type>,

            SpaceH::Interactions<BasicF, VelForce, ExtPosForce, ExtVelForce>,

            Regularitor<typename BasicF::type>
    >;

    template<
            typename BasicF      = SpaceH::NewtonianForce<DEFAULT_TYPE_CLASS>,
            typename VelForce    = void,
            typename ExtPosForce = void,
            typename ExtVelForce = void
    >
    using Basic = SpaceH::GARSystem<
            SpaceH::Particles<typename BasicF::type>,

            SpaceH::Interactions<BasicF, VelForce, ExtPosForce, ExtVelForce>,

            SpaceH::NoRegu<typename BasicF::type>
    >;

    /**  @brief Alias of template name, linking the particle system, integrator and ODE iterator*/
    template<
            typename ParticSys,
            template<typename, typename> class ODEiterator = SpaceH::BSIterator,
            template<typename> class Integrator  = SpaceH::symplectic2th
    >
    using Nbody = SpaceH::Solver<
            ParticSys,

            ODEiterator<ParticSys, Integrator<ParticSys> >
    >;
}
#endif
