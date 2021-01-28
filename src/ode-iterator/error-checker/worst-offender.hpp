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
 * @file worst-offender.hpp
 *
 * Header file.
 */
#pragma once

namespace space::ode_iterator {
    /*---------------------------------------------------------------------------*\
         Class WorstOffender Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam T
     */
    template <typename TypeSystem>
    class WorstOffender {
       public:
        // Type member
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(WorstOffender, default, default, default, default, default);

        WorstOffender(Scalar atol, Scalar rtol);

        SPACEHUB_READ_ACCESSOR(auto, atol, atol_);

        SPACEHUB_READ_ACCESSOR(auto, rtol, rtol_);

        void set_atol(Scalar);

        void set_rtol(Scalar);

        template <typename Array>
        auto error(Array const &y0, Array const &y1, Array const &y1_prime) -> typename Array::value_type;

        template <typename Array>
        auto error(Array const &y0, Array const &diff) -> typename Array::value_type;

       private:
        Scalar atol_{1e-13};

        Scalar rtol_{1e-13};
    };

    /*---------------------------------------------------------------------------*\
         Class WorstOffender Implementation
    \*---------------------------------------------------------------------------*/

    template <typename TypeSystem>
    void WorstOffender<TypeSystem>::set_atol(Scalar error) {
        atol_ = error;
    }

    template <typename TypeSystem>
    void WorstOffender<TypeSystem>::set_rtol(Scalar error) {
        rtol_ = error;
    }

    template <typename TypeSystem>
    template <typename Array>
    auto WorstOffender<TypeSystem>::error(const Array &y0, const Array &y1, const Array &y1_prime) ->
        typename Array::value_type {
        size_t const size = y0.size();
        Scalar max_err = 0;
        if constexpr (std::is_same_v<typename Array::value_type, Scalar>) {
            for (size_t i = 0; i < size; ++i) {
                Scalar scale = std::max(fabs(y1[i]), fabs(y0[i]));
                // Scalar scale = fabs(y0[i]) + fabs((y1[i] - y0[i]));
                max_err = math::max(max_err, fabs(y1_prime[i] - y1[i]) / (atol_ + scale * rtol_));
            }
        } else if constexpr (std::is_same_v<typename Array::value_type, Vec3<Scalar>>) {
            for (size_t i = 0; i < size; ++i) {
                Scalar scale = std::max(max_abs(y1[i]), max_abs(y0[i]));
                // Scalar scale = fabs(y0[i]) + fabs((y1[i] - y0[i]));
                max_err = math::max(max_err, max_abs(y1_prime[i] - y1[i]) / (atol_ + scale * rtol_));
            }
        } else {
            spacehub_abort("Unsupported array type!");
        }
        return max_err;
    }

    template <typename TypeSystem>
    template <typename Array>
    auto WorstOffender<TypeSystem>::error(const Array &y0, const Array &diff) -> typename Array::value_type {
        size_t const size = y0.size();
        Scalar max_err = 0;
        if constexpr (std::is_same_v<typename Array::value_type, Scalar>) {
            for (size_t i = 0; i < size; ++i) {
                Scalar scale = fabs(y0[i]);
                // Scalar scale = fabs(y0[i]) + fabs((y1[i] - y0[i]));
                max_err = math::max(max_err, fabs(diff[i]) / (atol_ + scale * rtol_));
            }
        } else if constexpr (std::is_same_v<typename Array::value_type, Vec3<Scalar>>) {
            for (size_t i = 0; i < size; ++i) {
                Scalar scale = max_abs(y0[i]);
                // Scalar scale = fabs(y0[i]) + fabs((y1[i] - y0[i]));
                max_err = math::max(max_err, max_abs(diff[i]) / (atol_ + scale * rtol_));
            }
        } else {
            spacehub_abort("Unsupported array type!");
        }
        return max_err;
    }
    template <typename TypeSystem>
    WorstOffender<TypeSystem>::WorstOffender(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}
}  // namespace space::ode_iterator
