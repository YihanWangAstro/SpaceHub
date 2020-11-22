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
 * @file IAS15-error.hpp
 *
 * Header file.
 */
#pragma once

namespace space::ode_iterator {
    /*---------------------------------------------------------------------------*\
         Class IAS15Error Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam T
     */
    template <typename T>
    class IAS15Error {
       public:
        // Type member

        using Scalar = T;

        using value_type = T;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(IAS15Error, default, default, default, default, default);

        IAS15Error(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

        SPACEHUB_READ_ACCESSOR(auto, atol, atol_);

        SPACEHUB_READ_ACCESSOR(auto, rtol, rtol_);

        void set_atol(Scalar);

        void set_rtol(Scalar);

        template <typename Array>
        auto error(Array const &scale, Array const &diff) -> typename Array::value_type;

        template <typename Array>
        auto error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type;

       private:
        Scalar atol_{1e-12};

        Scalar rtol_{1e-12};

        template <typename Array>
        auto one_dimension_error(Array const &scale, Array const &diff);

        template <typename Array>
        auto one_dimension_error(Array const &scale, Array const &y0, Array const &y1);
    };
    /*---------------------------------------------------------------------------*\
         Class IAS15Error Implementation
    \*---------------------------------------------------------------------------*/
    template <typename T>
    void IAS15Error<T>::set_atol(Scalar error) {
        atol_ = error;
    }

    template <typename T>
    void IAS15Error<T>::set_rtol(Scalar error) {
        rtol_ = error;
    }

    template <typename T>
    template <typename Array>
    auto IAS15Error<T>::error(const Array &scale, const Array &diff) -> typename Array::value_type {
        auto [max_diff, max_scale] = one_dimension_error(scale, diff);
        return max_diff / max_scale;
    }

    template <typename T>
    template <typename Array>
    auto IAS15Error<T>::error(const Array &scale, const Array &y0, const Array &y1) -> typename Array::value_type {
        auto [max_diff, max_scale] = one_dimension_error(scale, y0, y1);
        return max_diff / max_scale;
    }

    template <typename T>
    template <typename Array>
    auto IAS15Error<T>::one_dimension_error(const Array &scale, const Array &diff) {
        size_t const size = scale.size();
        Scalar max_diff = 0;
        Scalar max_scale = 0;
        for (size_t i = 0; i < size; ++i) {
            max_diff = std::max(max_diff, static_cast<Scalar>(fabs(diff[i])));
            max_scale = std::max(max_scale, static_cast<Scalar>(atol_ + fabs(scale[i]) * rtol_));
        }
        return std::make_tuple(max_diff, max_scale);
    }

    template <typename T>
    template <typename Array>
    auto IAS15Error<T>::one_dimension_error(const Array &scale, const Array &y0, const Array &y1) {
        size_t const size = scale.size();
        Scalar max_diff = 0;
        Scalar max_scale = 0;
        for (size_t i = 0; i < size; ++i) {
            max_diff = std::max(max_diff, static_cast<Scalar>(fabs(y0[i] - y1[i])));
            max_scale = std::max(max_scale, static_cast<Scalar>(atol_ + fabs(scale[i]) * rtol_));
        }
        return std::make_tuple(max_diff, max_scale);
    }
}  // namespace space::ode_iterator
