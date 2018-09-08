#ifndef SPACEHUB_H
#define SPACEHUB_H

#include "solver.h"
#include "particles.h"
#include "kahan_number.h"
#include "macros.h"

#include "tools/timer.h"
#include "callbacks.h"

#include "particle_system.h"
#include "particle_system/regularization.h"
#include "particle_system/regu_system.h"
#include "particle_system/GAR_system.h"

#include "particle_system/chain.h"
#include "particle_system/chain_particles.h"

#include "integrator/symplectic/symplectic2th.h"
#include "integrator/symplectic/symplectic4th.h"
#include "integrator/symplectic/symplectic6th.h"
#include "integrator/symplectic/symplectic8th.h"
#include "integrator/symplectic/symplectic10th.h"

#include "integrator/Gauss-Dadau.h"

#include "ode_iterator/BSIterator.h"
#include "ode_iterator/const_iterator.h"
#include "ode_iterator/IAS15.h"

#include "interaction/interaction.h"
#include "interaction/force.h"
#include "interaction/post_newtonian.h"

//#include "wrapper/singleRun.h"
namespace SpaceH {

    using DEFAULT_TYPE_CLASS = SpaceH::TypeClass<kahan<double>, SpaceH::DYNAMICAL>;

    template<
            typename BasicF      = SpaceH::NewtonianForce<DEFAULT_TYPE_CLASS>,
            typename VelForce    = void,
            typename ExtPosForce = void,
            typename ExtVelForce = void,
            template<typename> class Regularitor = SpaceH::LogH
    >
    using ARchain = SpaceH::GAR_system<
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
    using GAR = SpaceH::GAR_system<
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
    using Basic = SpaceH::GAR_system<
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
