/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @dir src @brief Root
 * @folder{vector}
 * @folder{tools}
 * @folder{stellar}
 * @folder{particles}
 * @folder{particle-system}
 * @folder{orbits}
 * @folder{ode-iterator}
 *   @folder{ode-iterator/error-checker}
 *   @folder{ode-iterator/step-controller}
 * @folder{multi-thread}
 * @folder{lazy-evaluation}
 * @folder{interaction}
 * @folder{integrator}
 *   @folder{integrator/symplectic}
 * @folder{args-callback}
 * @folder{scattering}
 *
 * @file spaceHub.hpp
 *
 * Giant header file.
 */

#pragma once

#ifdef MPFR_VERSION_MAJOR
#include "mpfr.hpp"
#endif

#include "args-callback/callbacks.hpp"
#include "integrator/Gauss-Radau.hpp"
#include "integrator/symplectic/symplectic-integrator.hpp"
#include "interaction/newtonian.hpp"
#include "interaction/post-newtonian.hpp"
#include "kahan-number.hpp"
#include "macros.hpp"
#include "multi-thread/multi-thread.hpp"
#include "ode-iterator/Bulirsch-Stoer.hpp"
#include "ode-iterator/IAS15.hpp"
#include "ode-iterator/const-iterator.hpp"
#include "ode-iterator/error-checker/RMS.hpp"
#include "ode-iterator/error-checker/max-raio-error.hpp"
#include "ode-iterator/error-checker/worst-offender.hpp"
#include "ode-iterator/sequent-iterator.hpp"
#include "ode-iterator/step-controller/PID-controller.hpp"
#include "ode-iterator/step-controller/const-controller.hpp"
#include "orbits/orbits.hpp"
#include "orbits/particle-manip.hpp"
#include "particle-system/archain.hpp"
#include "particle-system/base-system.hpp"
#include "particle-system/chain-system.hpp"
#include "particle-system/regu-system.hpp"
#include "particles/finite-size.hpp"
#include "particles/point-particles.hpp"
#include "scattering/cross-section.hpp"
#include "scattering/hierarchical.hpp"
#include "simulator.hpp"
#include "stellar/stellar.hpp"
#include "tools/auto-name.hpp"
#include "tools/config-reader.hpp"
#include "tools/timer.hpp"
#include "type-class.hpp"

/**
 * @namespace space
 * Documentation for space
 */
namespace space {

#define USING_NAMESPACE_SPACEHUB_ALL     \
    using namespace space;               \
    using namespace space::calc;         \
    using namespace space::tools;        \
    using namespace space::ode_iterator; \
    using namespace space::integrator;   \
    using namespace space::orbit;        \
    using namespace space::unit;         \
    using namespace space::particle_set; \
    using namespace space::random;       \
    using namespace space::callback;     \
    using namespace space::particle_system

    using DefaultTypes = Types<double, Vec3>;

    using DefaultForce = force::Interactions<space::force::NewtonianGrav>;

    template <typename T>
    using DefaultParticles = particle_set::PointParticles<T>;
    namespace methods {
        namespace details {
            using namespace ode_iterator;
            using namespace integrator;
            using normal_type = Types<double, Vec3>;
            using precise_type = Types<double_k, Vec3>;
#ifdef MPFR_VERSION_MAJOR
            using any_bits_type = Types<mpfr::mpreal, Vec3>;  // lazy vec3 will crash due to mpreal implementation.
            using precise_any_bits_type = Types<mpreal_k, Vec3>;
#endif
            using rms_err = ode_iterator::RMS<normal_type>;
            using worst_offender_err = ode_iterator::WorstOffender<normal_type>;
            using adaptive_step_ctrl = PIDController<normal_type>;
            using const_step_ctrl = ConstStepController<normal_type>;

            using const_sym2 = ConstOdeIterator<Symplectic2nd<normal_type>>;
            using const_sym4 = ConstOdeIterator<Symplectic4th<normal_type>>;
            using const_sym6 = ConstOdeIterator<Symplectic6th<normal_type>>;
            using const_sym8 = ConstOdeIterator<Symplectic8th<normal_type>>;
            using const_sym10 = ConstOdeIterator<Symplectic10th<normal_type>>;
            using const_Radau = ConstOdeIterator<GaussRadau<normal_type>>;

            using const_sym2_plus = ConstOdeIterator<Symplectic2nd<precise_type>>;
            using const_sym4_plus = ConstOdeIterator<Symplectic4th<precise_type>>;
            using const_sym6_plus = ConstOdeIterator<Symplectic6th<precise_type>>;
            using const_sym8_plus = ConstOdeIterator<Symplectic8th<precise_type>>;
            using const_sym10_plus = ConstOdeIterator<Symplectic10th<precise_type>>;
            using const_Radau_plus = ConstOdeIterator<GaussRadau<precise_type>>;

            using BS = BulirschStoer<LeapFrogDKD<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym2 = SequentOdeIterator<Symplectic2nd<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym4 = SequentOdeIterator<Symplectic4th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym6 = SequentOdeIterator<Symplectic6th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym8 = SequentOdeIterator<Symplectic8th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym10 = SequentOdeIterator<Symplectic10th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using Radau = IAS15<GaussRadau<normal_type>, MaxRatioError<normal_type>, adaptive_step_ctrl>;

            using BS_plus = BulirschStoer<LeapFrogDKD<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym2_plus = SequentOdeIterator<Symplectic2nd<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym4_plus = SequentOdeIterator<Symplectic4th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym6_plus = SequentOdeIterator<Symplectic6th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym8_plus = SequentOdeIterator<Symplectic8th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym10_plus = SequentOdeIterator<Symplectic10th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using Radau_plus = IAS15<GaussRadau<precise_type>, MaxRatioError<normal_type>, adaptive_step_ctrl>;
#ifdef MPFR_VERSION_MAJOR
            using ABits = BulirschStoer<LeapFrogDKD<any_bits_type>, ode_iterator::WorstOffender<any_bits_type>,
                                        PIDController<any_bits_type>, 32>;
            using ABits_plus =
                BulirschStoer<LeapFrogDKD<precise_any_bits_type>, ode_iterator::WorstOffender<precise_any_bits_type>,
                              PIDController<precise_any_bits_type>, 32>;
#endif
        };  // namespace details

#define DEFINE_ADAPTIVE_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                    \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using NAME = Simulator<particle_system::SYSTEM<particle<details::normal_type>, interactions>, details::ITER>; \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using NAME##_Plus =                                                                                           \
        Simulator<particle_system::SYSTEM<particle<details::precise_type>, interactions>, details::ITER##_plus>;

#define DEFINE_CONST_STEP_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME =                                                                                          \
        Simulator<particle_system::SYSTEM<particle<details::normal_type>, interactions>, details::const_##ITER>;  \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_Plus = Simulator<particle_system::SYSTEM<particle<details::precise_type>, interactions>, \
                                          details::const_##ITER##_plus>;

#define DEFINE_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                             \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using NAME = Simulator<particle_system::SYSTEM<particle<details::normal_type>, interactions>, details::ITER>; \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using NAME##_Plus =                                                                                           \
        Simulator<particle_system::SYSTEM<particle<details::precise_type>, interactions>, details::ITER##_plus>;  \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME =                                                                                          \
        Simulator<particle_system::SYSTEM<particle<details::normal_type>, interactions>, details::const_##ITER>;  \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_Plus = Simulator<particle_system::SYSTEM<particle<details::precise_type>, interactions>, \
                                          details::const_##ITER##_plus>;

#define DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(NAME, SYSTEM, ITER)                                                    \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>       \
    using NAME = Simulator<particle_system::SYSTEM<particle<details::any_bits_type>, interactions>, details::ITER>; \
                                                                                                                    \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>       \
    using NAME##_Plus = Simulator<particle_system::SYSTEM<particle<details::precise_any_bits_type>, interactions>,  \
                                  details::ITER##_plus>;

        /*template <typename force = DefaultForce, template <typename> typename particle = DefaultParticles>
        using AR_ABITS = Simulator<particle_system::RegularizedSystem<particle<details::any_bits_type>, force>,
                                   details::ABits>;

        template <typename force = DefaultForce, template <typename> typename particle = DefaultParticles>
        using ABITS =
            Simulator<particle_system::SimpleSystem<particle<details::any_bits_type>, force>, details::ABits>;*/

        DEFINE_ADAPTIVE_INTEGRATION_METHOD(BS, SimpleSystem, BS)

        DEFINE_ADAPTIVE_INTEGRATION_METHOD(AR_BS, RegularizedSystem, BS)

        DEFINE_ADAPTIVE_INTEGRATION_METHOD(Chain_BS, ChainSystem, BS)

        DEFINE_ADAPTIVE_INTEGRATION_METHOD(AR_Chain, ARchainSystem, BS)
#ifdef MPFR_VERSION_MAJOR
        DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(ABITS, SimpleSystem, ABits)

        DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(AR_ABITS, RegularizedSystem, ABits)
#endif
        DEFINE_INTEGRATION_METHOD(Sym2, SimpleSystem, sym2)

        DEFINE_INTEGRATION_METHOD(AR_Sym2, RegularizedSystem, sym2)

        DEFINE_INTEGRATION_METHOD(Chain_Sym2, ChainSystem, sym2)

        DEFINE_INTEGRATION_METHOD(AR_Sym2_Chain, ARchainSystem, sym2)

        DEFINE_INTEGRATION_METHOD(Sym4, SimpleSystem, sym4)

        DEFINE_INTEGRATION_METHOD(AR_Sym4, RegularizedSystem, sym4)

        DEFINE_INTEGRATION_METHOD(Chain_Sym4, ChainSystem, sym4)

        DEFINE_INTEGRATION_METHOD(AR_Sym4_Chain, ARchainSystem, sym4)

        DEFINE_INTEGRATION_METHOD(Sym6, SimpleSystem, sym6)

        DEFINE_INTEGRATION_METHOD(AR_Sym6, RegularizedSystem, sym6)

        DEFINE_INTEGRATION_METHOD(Chain_Sym6, ChainSystem, sym6)

        DEFINE_INTEGRATION_METHOD(AR_Sym6_Chain, ARchainSystem, sym6)

        DEFINE_INTEGRATION_METHOD(Sym8, SimpleSystem, sym8)

        DEFINE_INTEGRATION_METHOD(AR_Sym8, RegularizedSystem, sym8)

        DEFINE_INTEGRATION_METHOD(Chain_Sym8, ChainSystem, sym8)

        DEFINE_INTEGRATION_METHOD(AR_Sym8_Chain, ARchainSystem, sym8)

        DEFINE_INTEGRATION_METHOD(Sym10, SimpleSystem, sym10)

        DEFINE_INTEGRATION_METHOD(AR_Sym10, RegularizedSystem, sym10)

        DEFINE_INTEGRATION_METHOD(Chain_Sym10, ChainSystem, sym10)

        DEFINE_INTEGRATION_METHOD(AR_Sym10_Chain, ARchainSystem, sym10)

        DEFINE_INTEGRATION_METHOD(Radau, SimpleSystem, Radau)

        DEFINE_INTEGRATION_METHOD(AR_Radau, RegularizedSystem, Radau)

        DEFINE_INTEGRATION_METHOD(Chain_Radau, ChainSystem, Radau)

        DEFINE_INTEGRATION_METHOD(AR_Radau_Chain, ARchainSystem, Radau)
    }  // namespace methods

    using DefaultMethod = methods::AR_Chain_Plus<>;

    template <typename T>
    inline constexpr void set_mpreal_bits_from_rtol(T rtol) {
#ifdef MPFR_VERSION_MAJOR
        mpfr::mpreal::set_default_prec(size_t(mpfr::fabs(mpfr::LOG10(rtol))) * 4 + 32);
#endif
    }
}  // namespace space
