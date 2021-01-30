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
    template <typename TypeSystem>
    class IAS15Error {
       public:
        // Type member

        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(IAS15Error, default, default, default, default, default);

        IAS15Error(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

        SPACEHUB_READ_ACCESSOR(auto, atol, atol_);

        SPACEHUB_READ_ACCESSOR(auto, rtol, rtol_);

        void set_atol(Scalar);

        void set_rtol(Scalar);

        template <typename Array1, typename Array2>
        Scalar error(Array1 const &scale, Array2 const &diff);

        template <typename Array1, typename Array2, typename Array3>
        Scalar error(Array1 const &scale, Array2 const &y1, Array3 const &y1_prime);

       private:
        Scalar atol_{1e-12};

        Scalar rtol_{1e-12};
    };
    /*---------------------------------------------------------------------------*\
         Class IAS15Error Implementation
    \*---------------------------------------------------------------------------*/
    template <typename TypeSystem>
    void IAS15Error<TypeSystem>::set_atol(Scalar error) {
        atol_ = error;
    }

    template <typename TypeSystem>
    void IAS15Error<TypeSystem>::set_rtol(Scalar error) {
        rtol_ = error;
    }

    template <typename TypeSystem>
    template <typename Array1, typename Array2>
    auto IAS15Error<TypeSystem>::error(const Array1 &scale, const Array2 &diff) -> Scalar {
        size_t const size = scale.size();
        Scalar max_diff = 0;
        Scalar max_scale = 0;
        if constexpr (std::is_same_v<typename Array1::value_type, Scalar>) {
            for (size_t i = 0; i < size; ++i) {
                max_diff = std::max(max_diff, static_cast<Scalar>(fabs(diff[i])));
                max_scale = std::max(max_scale, static_cast<Scalar>(atol_ + fabs(scale[i]) * rtol_));
            }
        } else if constexpr (std::is_same_v<typename Array1::value_type, Vec3<Scalar>>) {
            for (size_t i = 0; i < size; ++i) {
                max_diff = std::max(max_diff, static_cast<Scalar>(max_abs(diff[i])));
                max_scale = std::max(max_scale, static_cast<Scalar>(atol_ + max_abs(scale[i]) * rtol_));
            }
        } else {
            spacehub_abort("Unsupported array type!");
        }
        return max_diff / max_scale;
    }

    template <typename TypeSystem>
    template <typename Array1, typename Array2, typename Array3>
    auto IAS15Error<TypeSystem>::error(const Array1 &scale, const Array2 &y1, const Array3 &y1_prime) -> Scalar {
        size_t const size = scale.size();
        Scalar max_diff = 0;
        Scalar max_scale = 0;
        if constexpr (std::is_same_v<typename Array1::value_type, Scalar>) {
            for (size_t i = 0; i < size; ++i) {
                max_diff = std::max(max_diff, static_cast<Scalar>(fabs(y1_prime[i] - y1[i])));
                max_scale = std::max(max_scale, static_cast<Scalar>(atol_ + fabs(scale[i]) * rtol_));
            }
        } else if constexpr (std::is_same_v<typename Array1::value_type, Vec3<Scalar>>) {
            for (size_t i = 0; i < size; ++i) {
                max_diff = std::max(max_diff, static_cast<Scalar>(max_abs(y1_prime[i] - y1[i])));
                max_scale = std::max(max_scale, static_cast<Scalar>(atol_ + max_abs(scale[i]) * rtol_));
            }
        } else {
            spacehub_abort("Unsupported array type!");
        }
        return max_diff / max_scale;
    }
}  // namespace space::ode_iterator
