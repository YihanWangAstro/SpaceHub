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
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file simulator.hpp
 *
 * Header file.
 */
#pragma once

#include <functional>

#include "IO.hpp"
#include "core-computation.hpp"
#include "dev-tools.hpp"

namespace space {

    /*---------------------------------------------------------------------------*\
        Class RunArgs Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Run arguments that is used to set all arguments needed by Simulator.
     *
     * The pre_operations(), post_operations(), stop_operations() and check_stops() will be called in a simulation
     * in this way.
     *
     * @code{.cpp}
     *  ...
     *  auto sys; // particle system that is going to be evolved.
     *  for( ;check_stops(sys);){
     *      pre_operations(sys);
     *      ...//step integration.
     *      post_operations(sys);
     *  }
     *  stop_operations(sys);
     *  ...
     * @endcode
     *
     * @tparam ParticleSys Any implementation of concept ParticleSystem.
     */
    template <typename ParticleSys>
    class RunArgs {
       public:
        // type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);

        /**
         * Callback function type for pre-operation, pos-operation and stop operation.
         */
        using Callback = std::function<void(ParticleSys &, Scalar)>;

        /**
         * Callback function type for stop condition.
         */
        using StopCall = std::function<bool(ParticleSys &, Scalar)>;

        // public members

        /**
         * Initial step size for the integration of the Simulator.
         */
        Scalar step_size{0};

        /**
         * The time duration for the integration of the Simulator.
         */
        Scalar end_time{0};

        /**
         * The absolute error tolerance.
         */
        Scalar atol{1e-14};

        /**
         * The relative error tolerance.
         */
        Scalar rtol{1e-13};

        // public methods
        /**
         * Call the all registered pre-operation functions by sequence.
         *
         * @param[in,out] particle_sys The particle system that is going to be operated.
         * @param[in] step_size The current step size.
         */
        void pre_operations(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Call the all registered pos-operation functions by sequence.
         *
         * @param[in,out] particle_sys The particle system that is going to be operated.
         * @param[in] step_size The current step size.
         */
        void post_operations(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Call the all registered stop-operation functions by sequence.
         *
         * @param[in,out] particle_sys The particle system that is going to be operated.
         * @param[in] step_size The current step size.
         */
        void stop_operations(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Check the all registered stop condition functions by sequence. If any of them is satisfied, return `true`.
         *
         * @param[in] particle_sys The particle system that is going to be checked
         * @param[in] step_size The current step size.
         */
        bool check_stops(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to pre-step-operations.
         *
         * @tparam Func Callable type that is conversional to member type Callback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments.If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_pre_step_operation(Func func, Args &&... args);

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to post-step-operations.
         *
         * @tparam Func Callable type that is conversional to member type Callback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments. If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_post_step_operation(Func func, Args &&... args);

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to stop-point-operations.
         *
         * @tparam Func Callable type that is conversional to member type Callback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments. If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_stop_point_operation(Func func, Args &&... args);

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to stop conditions.
         *
         * @tparam Func Callable type that is conversional to member type Stopback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments. If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_stop_condition(Func func, Args &&... args);

        /**
         * Add the duration time of the integration as a stop condition.
         * @tparam T Floating point like scalar.
         * @param[in] end Duration time of the integration.
         */
        template <typename T>
        void add_stop_condition(T end);

        /**
         * Check if the integration duration time is set.
         * @return boolean
         */
        [[nodiscard]] bool is_end_time_set() const { return is_end_time_set_; }

        /**
         * Check if any of the stop condition(except the duration time) is set.
         * @return boolean
         */
        [[nodiscard]] bool is_stop_condition_set() const { return stop_cond_.size() > 0; }

       private:
        // private members
        std::vector<Callback> pre_opts_;

        std::vector<Callback> post_opts_;

        std::vector<Callback> stop_opts_;

        std::vector<StopCall> stop_cond_;

        bool is_end_time_set_{false};
    };

    /*---------------------------------------------------------------------------*\
        Class Simulator Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Wrapper to integrate the particle system and ode-iterator to perform the simulation.
     *
     * @tparam ParticleSys Any implementation of concept `particle_system::ParticleSystem`.
     * @tparam OdeIterator Any implementation of concept `ode_iterator::OdeIterator`.
     */
    template <typename ParticleSys, typename OdeIterator>
    class Simulator {
       public:
        // Type member
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);

        /**
         * Run arguments that is used to set all arguments needed by Simulator.
         */
        using RunArgs = space::RunArgs<ParticleSys>;

        /**
         * Particle type that is used to create the initial conditions to initialize the Simulator.
         */
        using Particle = typename ParticleSys::Particle;

        SPACEHUB_READ_ACCESSOR(auto, particles, particles_.particles());

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(Simulator, delete, default, default, default, default);

        /**
         * Initialize the Simulator with an iterable Particle Container.
         * @tparam STL Iterable Particle Container.
         * @param[in] time Initial time of the particle system.
         * @param[in] particles_set Particle container.
         */
        template <typename STL>
        Simulator(Scalar time, STL const &particles_set);

        /**
         * Initialize the Simulator with an given particles.
         * @tparam T Any particle type that has the same interfaces of type member `Particle`.
         * @param[in] time Initial time of the particle system.
         * @param[in] particle Initial particles.
         */
        template <typename... T>
        explicit Simulator(Scalar time, T const &... particle);

        // Public methods
        /**
         * Run the simulation with given arguments.
         * @param[in] run_args Run arguments.
         */
        void run(RunArgs const &run_args);

        virtual ~Simulator() = default;

       private:
        // Private methods
        inline void advance_one_step();

        // Private members
        /** @brief Macro step size for ODE iterator*/
        Scalar step_size_{0.0};

        /** @brief Particle system*/
        ParticleSys particles_;

        /** @brief ODE Iterator*/
        OdeIterator iterator_;

        CREATE_METHOD_CHECK(set_atol);

        CREATE_METHOD_CHECK(set_rtol);
    };

    /*---------------------------------------------------------------------------*\
        Class RunArgs Implementation
    \*---------------------------------------------------------------------------*/
    template <typename ParticleSys>
    void RunArgs<ParticleSys>::pre_operations(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &opt : pre_opts_) {
            opt(particle_system, step_size);
        }
    }

    template <typename ParticleSys>
    void RunArgs<ParticleSys>::post_operations(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &opt : post_opts_) {
            opt(particle_system, step_size);
        }
    }

    template <typename ParticleSys>
    void RunArgs<ParticleSys>::stop_operations(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &opt : stop_opts_) {
            opt(particle_system, step_size);
        }
    }

    template <typename ParticleSys>
    bool RunArgs<ParticleSys>::check_stops(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &check : stop_cond_) {
            if (check(particle_system, step_size)) return true;
        }
        return false;
    }

    template <typename ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_pre_step_operation(Func func, Args &&... args) {
        pre_opts_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::placeholders::_2,
                                         std::forward<Args>(args)...));
    }

    template <typename ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_post_step_operation(Func func, Args &&... args) {
        post_opts_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::placeholders::_2,
                                          std::forward<Args>(args)...));
    }

    template <typename ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_stop_point_operation(Func func, Args &&... args) {
        stop_opts_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::placeholders::_2,
                                          std::forward<Args>(args)...));
    }

    template <typename ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_stop_condition(Func func, Args &&... args) {
        stop_cond_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::placeholders::_2,
                                          std::forward<Args>(args)...));
    }

    template <typename ParticleSys>
    template <typename T>
    void RunArgs<ParticleSys>::add_stop_condition(T end) {
        end_time = end;
        is_end_time_set_ = true;
    }

    /*---------------------------------------------------------------------------*\
        Class Simulator Implememtation
    \*---------------------------------------------------------------------------*/
    template <typename ParticleSys, typename OdeIterator>
    template <typename STL>
    Simulator<ParticleSys, OdeIterator>::Simulator(Scalar time, const STL &particle_set)
        : particles_(time, particle_set) {
        static_assert(is_ranges_v<STL>, "Only STL-like container can be used");
    }

    template <typename ParticleSys, typename OdeIterator>
    template <typename... T>
    Simulator<ParticleSys, OdeIterator>::Simulator(Scalar time, T const &... particle)
        : Simulator(time, std::initializer_list<Particle>{particle...}) {
        static_assert(calc::all(std::is_same_v<T, Particle>...), "Wrong particles type!");
    }

    template <typename ParticleSys, typename OdeIterator>
    void Simulator<ParticleSys, OdeIterator>::run(RunArgs const &run_args) {
        if (!run_args.is_stop_condition_set() && !run_args.is_end_time_set()) {
            space::spacehub_abort("Use 'add_stop_condition' to set stop condition.");
        }

        step_size_ = run_args.step_size;

        if (step_size_ == 0.0) {
            step_size_ = 0.01 * calc::calc_step_scale(particles_.particles());
        }

        Scalar end_time = space::unit::T_hubble;

        if (run_args.is_end_time_set()) {
            end_time = run_args.end_time;
        }

        if (particles_.particles().time() >= end_time) {
            space::print(std::cout, "Warning: The stop time is '<=' to the start time!");
        }

        if constexpr (HAS_METHOD(OdeIterator, set_atol, Scalar)) {
            iterator_.set_atol(run_args.atol);
        }

        if constexpr (HAS_METHOD(OdeIterator, set_rtol, Scalar)) {
            iterator_.set_rtol(run_args.rtol);
        }

        for (; particles_.particles().time() < end_time && !run_args.check_stops(particles_, step_size_);) {
            run_args.pre_operations(particles_, step_size_);
            advance_one_step();
            run_args.post_operations(particles_, step_size_);
        }
        run_args.stop_operations(particles_, step_size_);
    }

    template <typename ParticleSys, typename OdeIterator>
    inline void Simulator<ParticleSys, OdeIterator>::advance_one_step() {
        particles_.pre_iter_process();
        step_size_ = iterator_.iterate(particles_, step_size_);
        particles_.post_iter_process();
    }

}  // namespace space
