#ifndef SPACEX_H
#define SPACEX_H
#include "dynamicSystem.h"
#include "dynamicState.h"
#include "integrator/symplectic/symplectic2th.h"
#include "ODEiterator/BSIterator.h"
#include "ODEiterator/constIterator.h"
#include "particleSystem/reguSystem.h"
#include "particleSystem/ARchain.h"
#include "particleSystem/GAR.h"
#include "particleSystem/regularState.h"
#include "particleSystem/regularization.h"
#include "interaction/interaction.h"

template<
size_t                            N,
template<typename> class          Regularitor = logH,
typename                          Scalar      = double,
template<typename, size_t> class  EvolvedData = reguDynamics,
template<typename> class          Interaction = Newtonian
>
using NewtonianSystem = reguSystem<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;

template<
size_t                            N,
template<typename> class          Interaction,
template<typename> class          Regularitor = logH,
typename                          Scalar      = double,
template<typename, size_t> class  EvolvedData = GAR
>
using VelDepSystem = reguSystem<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;


template<
size_t                            N,
template<typename> class          Regularitor = logH,
typename                          Scalar      = double ,
template<typename, size_t> class  EvolvedData = GAR,
template<typename> class          Interaction = Newtonian
>
using AR_chain = reguSystem<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;

template<
size_t                            N,
template<typename> class          Interaction,
template<typename> class          Regularitor = logH,
typename                          Scalar      = double,
template<typename, size_t> class  EvolvedData = GAR
>
using VelARchain = ARchain<Interaction<Scalar>, EvolvedData<Scalar, N>, Regularitor<EvolvedData<Scalar, N>>>;

#endif
