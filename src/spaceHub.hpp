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
#include "args-callback/collision.hpp"
#include "integrator/Gauss-Radau.hpp"
#include "integrator/symplectic/symplectic-integrator.hpp"
// #include "interaction/disk-capture-star.hpp"
// #include "interaction/disk-capture.hpp"
#include "interaction/magneto-disk.hpp"
#include "interaction/newtonian.hpp"
#include "interaction/post-newtonian.hpp"
#include "interaction/tidal.hpp"
#include "kahan-number.hpp"
#include "macros.hpp"
#include "multi-thread/multi-thread.hpp"
#include "ode-iterator/Bulirsch-Stoer.hpp"
#include "ode-iterator/IAS15.hpp"
#include "ode-iterator/const-iterator.hpp"
#include "ode-iterator/error-checker/RMS.hpp"
#include "ode-iterator/error-checker/max-ratio-error.hpp"
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
#include "particles/drag-particles.hpp"
#include "particles/finite-size.hpp"
#include "particles/point-particles.hpp"
#include "particles/tide-particles.hpp"
#include "scattering/cross-section.hpp"
#include "scattering/hierarchical.hpp"
#include "simulator.hpp"
#include "stellar/stellar.hpp"
#include "tools/auto-name.hpp"
#include "tools/config-reader.hpp"
#include "tools/timer.hpp"
#include "type-class.hpp"

/**
 * @namespace hub
 * Documentation for hub
 */
namespace hub {

#define USING_NAMESPACE_SPACEHUB_ALL \
    using namespace hub;             \
    using namespace hub::calc;       \
    using namespace hub::tools;      \
    using namespace hub::ode;        \
    using namespace hub::integrator; \
    using namespace hub::orbit;      \
    using namespace hub::unit;       \
    using namespace hub::particles;  \
    using namespace hub::random;     \
    using namespace hub::callback;   \
    using namespace hub::system;     \
    using namespace hub::force

    using DefaultTypes = Types<double, Vec3>;

    using DefaultForce = force::Interactions<hub::force::NewtonianGrav>;

    template <typename T>
    using DefaultParticles = particles::PointParticles<T>;
    namespace methods {
        namespace details {
            using namespace ode;
            using namespace integrator;
            using normal_type = Types<double, Vec3>;
            using extended_type = Types<long double, Vec3>;
            using precise_type = Types<double_k, Vec3>;
            using extended_precise_type = Types<long_double_k, Vec3>;
#ifdef MPFR_VERSION_MAJOR
            using any_bits_type = Types<mpfr::mpreal, Vec3>;  // lazy vec3 will crash due to mpreal implementation.
            using precise_any_bits_type = Types<mpreal_k, Vec3>;
#endif
            using rms_err = ode::RMS<normal_type>;
            using worst_offender_err = ode::WorstOffender<normal_type>;
            using adaptive_step_ctrl = PIDController<normal_type>;
            using const_step_ctrl = ConstStepController<normal_type>;

            using rms_err_ext = ode::RMS<extended_type>;
            using worst_offender_err_ext = ode::WorstOffender<extended_type>;
            using adaptive_step_ctrl_ext = PIDController<extended_type>;
            using const_step_ctrl_ext = ConstStepController<extended_type>;

            using const_sym2 = ConstOdeIterator<Symplectic2nd<normal_type>>;
            using const_sym4 = ConstOdeIterator<Symplectic4th<normal_type>>;
            using const_sym6 = ConstOdeIterator<Symplectic6th<normal_type>>;
            using const_sym8 = ConstOdeIterator<Symplectic8th<normal_type>>;
            using const_sym10 = ConstOdeIterator<Symplectic10th<normal_type>>;
            using const_Radau = ConstOdeIterator<GaussRadau<normal_type>>;

            using const_sym2_ext = ConstOdeIterator<Symplectic2nd<extended_type>>;
            using const_sym4_ext = ConstOdeIterator<Symplectic4th<extended_type>>;
            using const_sym6_ext = ConstOdeIterator<Symplectic6th<extended_type>>;
            using const_sym8_ext = ConstOdeIterator<Symplectic8th<extended_type>>;
            using const_sym10_ext = ConstOdeIterator<Symplectic10th<extended_type>>;
            using const_Radau_ext = ConstOdeIterator<GaussRadau<extended_type>>;

            using const_sym2_plus = ConstOdeIterator<Symplectic2nd<precise_type>>;
            using const_sym4_plus = ConstOdeIterator<Symplectic4th<precise_type>>;
            using const_sym6_plus = ConstOdeIterator<Symplectic6th<precise_type>>;
            using const_sym8_plus = ConstOdeIterator<Symplectic8th<precise_type>>;
            using const_sym10_plus = ConstOdeIterator<Symplectic10th<precise_type>>;
            using const_Radau_plus = ConstOdeIterator<GaussRadau<precise_type>>;

            using const_sym2_extplus = ConstOdeIterator<Symplectic2nd<extended_precise_type>>;
            using const_sym4_extplus = ConstOdeIterator<Symplectic4th<extended_precise_type>>;
            using const_sym6_extplus = ConstOdeIterator<Symplectic6th<extended_precise_type>>;
            using const_sym8_extplus = ConstOdeIterator<Symplectic8th<extended_precise_type>>;
            using const_sym10_extplus = ConstOdeIterator<Symplectic10th<extended_precise_type>>;
            using const_Radau_extplus = ConstOdeIterator<GaussRadau<extended_precise_type>>;

            using BS = BulirschStoer<LeapFrogDKD<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym2 = SequentOdeIterator<Symplectic2nd<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym4 = SequentOdeIterator<Symplectic4th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym6 = SequentOdeIterator<Symplectic6th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym8 = SequentOdeIterator<Symplectic8th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym10 = SequentOdeIterator<Symplectic10th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using Radau = IAS15<GaussRadau<normal_type>, MaxRatioError<normal_type>, adaptive_step_ctrl>;

            using BS_ext = BulirschStoer<LeapFrogDKD<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym2_ext =
                SequentOdeIterator<Symplectic2nd<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym4_ext =
                SequentOdeIterator<Symplectic4th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym6_ext =
                SequentOdeIterator<Symplectic6th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym8_ext =
                SequentOdeIterator<Symplectic8th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym10_ext =
                SequentOdeIterator<Symplectic10th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using Radau_ext = IAS15<GaussRadau<extended_type>, MaxRatioError<extended_type>, adaptive_step_ctrl_ext>;

            using BS_plus = BulirschStoer<LeapFrogDKD<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym2_plus = SequentOdeIterator<Symplectic2nd<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym4_plus = SequentOdeIterator<Symplectic4th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym6_plus = SequentOdeIterator<Symplectic6th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym8_plus = SequentOdeIterator<Symplectic8th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym10_plus = SequentOdeIterator<Symplectic10th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using Radau_plus = IAS15<GaussRadau<precise_type>, MaxRatioError<normal_type>, adaptive_step_ctrl>;

            using BS_extplus =
                BulirschStoer<LeapFrogDKD<extended_precise_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym2_extplus = SequentOdeIterator<Symplectic2nd<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym4_extplus = SequentOdeIterator<Symplectic4th<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym6_extplus = SequentOdeIterator<Symplectic6th<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym8_extplus = SequentOdeIterator<Symplectic8th<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym10_extplus = SequentOdeIterator<Symplectic10th<extended_precise_type>, worst_offender_err_ext,
                                                     adaptive_step_ctrl_ext>;
            using Radau_extplus =
                IAS15<GaussRadau<extended_precise_type>, MaxRatioError<extended_type>, adaptive_step_ctrl_ext>;
#ifdef MPFR_VERSION_MAJOR
            using ABits = BulirschStoer<LeapFrogDKD<any_bits_type>, ode::WorstOffender<any_bits_type>,
                                        PIDController<any_bits_type>, 32>;
            using ABits_plus =
                BulirschStoer<LeapFrogDKD<precise_any_bits_type>, ode::WorstOffender<precise_any_bits_type>,
                              PIDController<precise_any_bits_type>, 32>;
#endif
        };  // namespace details

#define DEFINE_ADAPTIVE_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                         \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME = Simulator<system::SYSTEM<particle<details::normal_type>, interactions>, details::ITER>;               \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Plus =                                                                                                \
        Simulator<system::SYSTEM<particle<details::precise_type>, interactions>, details::ITER##_plus>;                \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Ext = Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::ITER##_ext>; \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_ExtPlus =                                                                                             \
        Simulator<system::SYSTEM<particle<details::extended_precise_type>, interactions>, details::ITER##_extplus>;

#define DEFINE_CONST_STEP_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME =                                                                                          \
        Simulator<particle_system::SYSTEM<particle<details::normal_type>, interactions>, details::const_##ITER>;  \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_Plus = Simulator<particle_system::SYSTEM<particle<details::precise_type>, interactions>, \
                                          details::const_##ITER##_plus>;                                          \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_Ext = Simulator<particle_system::SYSTEM<particle<details::extended_type>, interactions>, \
                                         details::const_##ITER##_ext>;                                            \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_ExtPlus =                                                                                \
        Simulator<particle_system::SYSTEM<particle<details::extended_precise_type>, interactions>,                \
                  details::const_##ITER##_extplus>;

#define DEFINE_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME = Simulator<system::SYSTEM<particle<details::normal_type>, interactions>, details::ITER>;               \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Plus =                                                                                                \
        Simulator<system::SYSTEM<particle<details::precise_type>, interactions>, details::ITER##_plus>;                \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME =                                                                                               \
        Simulator<system::SYSTEM<particle<details::normal_type>, interactions>, details::const_##ITER>;                \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME##_Plus =                                                                                        \
        Simulator<system::SYSTEM<particle<details::precise_type>, interactions>, details::const_##ITER##_plus>;        \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Ext = Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::ITER##_ext>; \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_ExtPlus =                                                                                             \
        Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::ITER##_extplus>;            \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME##_Ext =                                                                                         \
        Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::const_##ITER##_ext>;        \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME##_ExtPlus = Simulator<system::SYSTEM<particle<details::extended_precise_type>, interactions>,   \
                                             details::const_##ITER##_extplus>;

#define DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(NAME, SYSTEM, ITER)                                              \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles> \
    using NAME = Simulator<system::SYSTEM<particle<details::any_bits_type>, interactions>, details::ITER>;    \
                                                                                                              \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles> \
    using NAME##_Plus =                                                                                       \
        Simulator<system::SYSTEM<particle<details::precise_any_bits_type>, interactions>, details::ITER##_plus>;

        /*template <typename force = DefaultForce, template <typename> typename particle = DefaultParticles>
        using AR_ABITS = Simulator<system::RegularizedSystem<particle<details::any_bits_type>, force>,
                                   details::ABits>;

        template <typename force = DefaultForce, template <typename> typename particle = DefaultParticles>
        using ABITS =
            Simulator<system::SimpleSystem<particle<details::any_bits_type>, force>, details::ABits>;*/

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

        template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>
        using DefaultMethod = methods::AR_Chain_Plus<interactions, particle>;
    }  // namespace methods

    template <typename T>
    inline constexpr void set_mpreal_bits_from_rtol(T rtol) {
#ifdef MPFR_VERSION_MAJOR
        mpfr::mpreal::set_default_prec(size_t(mpfr::fabs(mpfr::LOG10(rtol))) * 4 + 32);
#endif
    }
}  // namespace hub
