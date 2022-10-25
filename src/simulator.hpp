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
 * @file simulator.hpp
 *
 * Header file.
 */
#pragma once
#include <functional>

#include "IO.hpp"
#include "core-computation.hpp"
#include "dev-tools.hpp"
namespace hub {

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
    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
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
        Scalar atol{0.0};

        /**
         * The relative error tolerance.
         */
        Scalar rtol{1e-14};

        /** Added by Dipto
         * Relative time tolerance parameter to iterate to exact time for
         * regularized methods
         */

        Scalar time_rtol{1e-8};

        // public methods
        /**
         * Call the all registered operation functions by sequence.
         *
         * @param[in,out] particle_sys The particle system that is going to be operated.
         * @param[in] step_size The current step size.
         */
        void operations(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Call the all registered stop-operation functions by sequence.
         *
         * @param[in,out] particle_sys The particle system that is going to be operated.
         * @param[in] step_size The current step size.
         */
        void stop_operations(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Call the all registered start-operation functions by sequence.
         *
         * @param[in,out] particle_sys The particle system that is going to be operated.
         * @param[in] step_size The current step size.
         */
        void start_operations(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Check the all registered stop condition functions by sequence. If any of them is satisfied, return `true`.
         *
         * @param[in] particle_sys The particle system that is going to be checked
         * @param[in] step_size The current step size.
         */
        bool check_stops(ParticleSys &particle_sys, Scalar step_size) const;

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to start point and every post-step.
         *
         * @tparam Func Callable type that is convertible to member type Callback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments.If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_operation(Func func, Args &&...args);

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to stop-point-operations.
         *
         * @tparam Func Callable type that is convertible to member type Callback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments. If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_stop_point_operation(Func func, Args &&...args);

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to start-point-operations.
         *
         * @tparam Func Callable type that is convertible to member type Callback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments. If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_start_point_operation(Func func, Args &&...args);

        /**
         * Register a callable object(function pointer, functor, lambda,etc...) to stop conditions.
         *
         * @tparam Func Callable type that is convertible to member type Stopback.
         * @tparam Args Type of the binding arguments.
         * @param[in] func Callable object.
         * @param[in] args Binding arguments. If func accepts more than one arguments, you can bind the rest arguments
         * here.
         */
        template <typename Func, typename... Args>
        void add_stop_condition(Func func, Args &&...args);

        /**
         * Add the duration time of the integration as a stop condition.
         * @param[in] end Duration time of the integration.
         */
        void add_stop_condition(Scalar end);

        /**
         * Check if the integration duration time is set.
         * @return boolean
         */
        bool is_end_time_set() const { return is_end_time_set_; }

        /**
         * Check if any of the stop condition(except the duration time) is set.
         * @return boolean
         */
        bool is_stop_condition_set() const { return stop_cond_.size() > 0; }

       private:
        // private members
        std::vector<Callback> opts_;

        std::vector<Callback> stop_opts_;

        std::vector<Callback> start_opts_;

        std::vector<StopCall> stop_cond_;

        bool is_end_time_set_{false};

        CREATE_METHOD_CHECK(operation);
    };

    /*---------------------------------------------------------------------------*\
        Class Simulator Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Wrapper to integrate the particle system and ode-iterator to perform the simulation.
     *
     * @tparam ParticleSys Any implementation of concept `system::ParticleSystem`.
     * @tparam OdeIterator Any implementation of concept `ode::OdeIterator`.
     */
    template <typename ParticleSys, typename OdeIterator>
    class Simulator {
       public:
        // Type member
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);

        /**
         * Run arguments that is used to set all arguments needed by Simulator.
         */
        using RunArgs = hub::RunArgs<ParticleSys>;

        using ParticleSystem = ParticleSys;
 
        /**
         * Particle type that is used to create the initial conditions to initialize the Simulator.
         */
        using Particle = typename ParticleSys::Particle;

        SPACEHUB_READ_ACCESSOR(ParticleSys, particles, particles_);

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(Simulator, delete, default, default, default, default);

        /**
         * Initialize the Simulator with an iterable Particle Container.
         * @tparam STL Iterable Particle Container.
         * @param[in] time Initial time of the particle system.
         * @param[in] particles_set Particle container.
         */
        template <CONCEPT_PARTICLE_CONTAINER STL>
        Simulator(Scalar time, STL const &particles_set);

        /**
         * Initialize the Simulator with an given particles.
         * @tparam T Any particle type that has the same interfaces of type member `Particle`.
         * @param[in] time Initial time of the particle system.
         * @param[in] particle Initial particles.
         */
        template <typename... T>
        explicit Simulator(Scalar time, T const &...particle);

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
    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    void RunArgs<ParticleSys>::operations(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &opt : opts_) {
            opt(particle_system, step_size);
        }
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    void RunArgs<ParticleSys>::stop_operations(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &opt : stop_opts_) {
            opt(particle_system, step_size);
        }
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    void RunArgs<ParticleSys>::start_operations(ParticleSys &particle_system, Scalar step_size) const {
        for (auto const &opt : start_opts_) {
            opt(particle_system, step_size);
        }
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    bool RunArgs<ParticleSys>::check_stops(ParticleSys &particle_system, Scalar step_size) const {
        return std::any_of(stop_cond_.begin(), stop_cond_.end(),
                           [&](auto cond) { return cond(particle_system, step_size); });
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_operation(Func func, Args &&...args) {
        opts_.emplace_back(std::bind(func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...));
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_stop_point_operation(Func func, Args &&...args) {
        stop_opts_.emplace_back(
            std::bind(func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...));
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_start_point_operation(Func func, Args &&...args) {
        start_opts_.emplace_back(
            std::bind(func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...));
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    template <typename Func, typename... Args>
    void RunArgs<ParticleSys>::add_stop_condition(Func func, Args &&...args) {
        stop_cond_.emplace_back(
            std::bind(func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...));
    }

    template <CONCEPT_PARTICLE_SYSTEM ParticleSys>
    void RunArgs<ParticleSys>::add_stop_condition(Scalar end) {
        end_time = end;
        is_end_time_set_ = true;
    }

    /*---------------------------------------------------------------------------*\
        Class Simulator Implementation
    \*---------------------------------------------------------------------------*/
    template <typename ParticleSys, typename OdeIterator>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    Simulator<ParticleSys, OdeIterator>::Simulator(Scalar time, const STL &particle_set)
        : particles_(time, particle_set) {}

    template <typename ParticleSys, typename OdeIterator>
    template <typename... T>
    Simulator<ParticleSys, OdeIterator>::Simulator(Scalar time, T const &...particle)
        : Simulator(time, std::initializer_list<Particle>{particle...}) {
        static_assert(calc::all(std::is_same_v<T, Particle>...), "Wrong particles type!");
    }

    template <typename ParticleSys, typename OdeIterator>
    void Simulator<ParticleSys, OdeIterator>::run(RunArgs const &run_args) {
        if (!run_args.is_stop_condition_set() && !run_args.is_end_time_set()) {
            spacehub_abort("Use 'add_stop_condition' to set stop condition.");
        }

        step_size_ = run_args.step_size;
        auto time_rtol_ = run_args.time_rtol;
        if (step_size_ == 0.0) {
            step_size_ = 0.1 * calc::calc_step_scale(particles_) *
                         calc::calc_fall_free_time(particles_.mass(), particles_.pos());
        }

        Scalar end_time = hub::unit::T_hubble;

        if (run_args.is_end_time_set()) {
            end_time = run_args.end_time;
        }

        if (particles_.time() >= end_time) {
            hub::print(std::cout, "Warning: The stop time is '<=' to the start time!");
        }

        if (particles_.step_scale() == 0.0) {
            hub::print(std::cout, "regularization function === 0. Use other method instead\n");
            return;
        }

        if constexpr (HAS_METHOD(OdeIterator, set_atol, Scalar)) {
            iterator_.set_atol(run_args.atol);
        }

        if constexpr (HAS_METHOD(OdeIterator, set_rtol, Scalar)) {
            iterator_.set_rtol(run_args.rtol);
        }

        run_args.start_operations(particles_, step_size_);

        // Dipto's changes here

        for (; std::abs((particles_.time() - end_time) / end_time) > time_rtol_ &&
               !run_args.check_stops(particles_, step_size_);) {
            Scalar rest_step = (end_time - particles_.time()) * particles_.step_scale();

            if (std::abs(step_size_) <= std::abs(rest_step)) [[likely]] {
                run_args.operations(particles_, step_size_);
                advance_one_step();

            } else {
                step_size_ = rest_step;
                run_args.operations(particles_, step_size_);
                advance_one_step();
            }
        }
        run_args.operations(particles_, step_size_);
        run_args.stop_operations(particles_, step_size_);
    }

    template <typename ParticleSys, typename OdeIterator>
    inline void Simulator<ParticleSys, OdeIterator>::advance_one_step() {
        particles_.pre_iter_process();
        step_size_ = iterator_.iterate(particles_, step_size_);
        particles_.post_iter_process();
    }

}  // namespace hub
