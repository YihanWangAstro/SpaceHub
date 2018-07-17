#ifndef SPACEX_H
#define SPACEX_H
#include "dynamicSystem.h"
#include "particles.h"
#include "kahanNumber.h"

#include "particleSystem.h"
#include "particleSystem/reguSystem/regularization.h"
#include "particleSystem/reguSystem/reguSystem.h"
#include "particleSystem/reguSystem/regularState.h"

//#include "particleSystem/ARchain/ARchain.h"
#include "particleSystem/ARChain/chain.h"
#include "particleSystem/ARChain/dynamicChain.h"

#include "integrator/symplectic/symplectic2th.h"
#include "integrator/symplectic/symplectic4th.h"
#include "integrator/symplectic/symplectic6th.h"
#include "integrator/symplectic/symplectic8th.h"
#include "integrator/symplectic/symplectic10th.h"
#include "ODEiterator/BSIterator.h"
#include "ODEiterator/constIterator.h"

#include "interaction/interaction.h"
#include "interaction/forces.h"
#include "interaction/postNewtonian.h"
/*
template <
    size_t                            N,
    template<typename> class          Regularitor = logH,
    typename                          Scalar      = double,
    template<typename, size_t> class  EvolvedData = reguDynamics,
template<typename> class          Interaction = interact::Newtonian
    >
using NewtonianSystem = reguSystem<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;

template <
    size_t                            N,
    template<typename> class          Interaction,
    template<typename> class          Regularitor = logH,
    typename                          Scalar      = double,
    template<typename, size_t> class  EvolvedData = GAR
    >
using VelDepSystem = reguSystem<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;


template <
    size_t                            N,
    template<typename> class          Regularitor = logH,
    typename                          Scalar      = double,
    template<typename, size_t> class  EvolvedData = reguDynamics,
    template<typename> class          Interaction = interact::Newtonian
    >
using AR_chain = ARchain<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;

template <
    size_t                            N,
    template<typename> class          Interaction,
    template<typename> class          Regularitor = logH,
    typename                          Scalar      = double,
    template<typename, size_t> class  EvolvedData = GAR
    >
using VelARchain = ARchain<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;
*/
#endif
