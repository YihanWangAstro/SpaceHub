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
 * @file RMS.hpp
 *
 * Header file.
 */
#pragma once

namespace space::ode_iterator {
    /*---------------------------------------------------------------------------*\
         Class RMS Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam T
     */
    template <typename T>
    class RMS {
       public:
        // Type member

        using Scalar = T;

        using value_type = T;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(RMS, default, default, default, default, default);

        RMS(Scalar atol, Scalar rtol);

        SPACEHUB_READ_ACCESSOR(auto, atol, atol_);

        SPACEHUB_READ_ACCESSOR(auto, rtol, rtol_);

        void set_atol(Scalar);

        void set_rtol(Scalar);

        template <typename Array>
        auto error(Array const &y0, Array const &y1, Array const &y1_prime) -> typename Array::value_type;

       private:
        Scalar atol_{1e-13};

        Scalar rtol_{1e-13};
    };
    /*---------------------------------------------------------------------------*\
         Class RMS Implementation
    \*---------------------------------------------------------------------------*/
    template <typename T>
    void RMS<T>::set_atol(Scalar error) {
        atol_ = error;
    }

    template <typename T>
    void RMS<T>::set_rtol(Scalar error) {
        rtol_ = error;
    }

    template <typename T>
    template <typename Array>
    auto RMS<T>::error(const Array &y0, const Array &y1, const Array &y1_prime) -> typename Array::value_type {
        size_t const size = y0.size();
        Scalar error = 0;
        for (size_t i = 0; i < size; ++i) {
            Scalar scale = std::max(fabs(y0[i]), fabs(y1[i]));
            auto r = fabs(y1[i] - y1_prime[i]) / (atol_ + scale * rtol_);
            error += r * r;
        }
        return sqrt(error / size);
    }

    template <typename T>
    RMS<T>::RMS(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}
}  // namespace space::ode_iterator
