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
 * @file const-iterator.hpp
 *
 * Header file.
 */
#pragma once
#include "../spacehub-concepts.hpp"
namespace space::ode_iterator {

    /*---------------------------------------------------------------------------*\
          Class BisecOdeIterator Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam Integrator
     */
    template <typename Integrator, typename ErrEstimator, typename StepController>
    class BisecOdeIterator {
       public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Integrator);

        using State = ScalarArray;

        template <CONCEPT_PARTICLE_SYSTEM T>
        auto iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar;

       private:
        void check_variable_size();

        // Private members
        Integrator integrator_;

        ErrEstimator err_checker_;

        StepController step_controller_;

        State input_;

        State output_;

        State dual_steps_output_;

        size_t var_num_;
    };

    /*---------------------------------------------------------------------------*\
          Class BisecOdeIterator Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Integrator, typename ErrEstimator, typename StepController>
    template <CONCEPT_PARTICLE_SYSTEM T>
    auto BisecOdeIterator<Integrator, ErrEstimator, StepController>::iterate(T &particles,
                                                                             typename T::Scalar macro_step_size) ->
        typename T::Scalar {
        check_variable_size();
        particles.write_to_scalar_array(input_);

        size_t n_steps = 1;

        Scalar h = macro_step_size;

        Scalar error = 0;

        integrator_.integrate(particles, h);

        particles.write_to_scalar_array(output_);

        for (;;) {
            h /= 2;
            n_steps *= 2;
            particles.read_from_scalar_array(input_);
            for (size_t i = 0; i < n_steps; ++i) {
                integrator_.integrate(particles, h);
            }
            particles.write_to_scalar_array(dual_steps_output_);
            error = err_checker_.error(input_, output_, dual_steps_output_);
            if (error <= 1) {
                break;
            }
            std::swap(output_, dual_steps_output_);
        }
        particles.read_from_scalar_array(dual_steps_output_);
        return step_controller_.next_step_size(Integrator::order, macro_step_size, error);
    }
    template <typename Integrator, typename ErrEstimator, typename StepController>
    void BisecOdeIterator<Integrator, ErrEstimator, StepController>::check_variable_size() {
        var_num_ = input_.size();
        if (var_num_ > output_.size()) [[unlikely]] {
            output_.resize(var_num);
            dual_steps_output_.resize(var_num);
        }
    }
}  // namespace space::ode_iterator
