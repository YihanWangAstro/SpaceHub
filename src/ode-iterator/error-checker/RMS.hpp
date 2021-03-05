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
 * @file RMS.hpp
 *
 * Header file.
 */
#pragma once

namespace hub::ode {
    /*---------------------------------------------------------------------------*\
         Class RMS Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief
     *
     * @tparam TypeSystem
     */
    template <typename TypeSystem>
    class RMS {
       public:
        // Type member

        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(RMS, default, default, default, default, default);

        RMS(Scalar atol, Scalar rtol);

        SPACEHUB_READ_ACCESSOR(Scalar, atol, atol_);

        SPACEHUB_READ_ACCESSOR(Scalar, rtol, rtol_);

        void set_atol(Scalar);

        void set_rtol(Scalar);

        template <typename Array1, typename Array2, typename Array3>
        Scalar error(Array1 const &y0, Array2 const &y1, Array3 const &y1_prime);

       private:
        Scalar atol_{1e-13};

        Scalar rtol_{1e-13};
    };
    /*---------------------------------------------------------------------------*\
         Class RMS Implementation
    \*---------------------------------------------------------------------------*/
    template <typename TypeSystem>
    void RMS<TypeSystem>::set_atol(Scalar error) {
        atol_ = error;
    }

    template <typename TypeSystem>
    void RMS<TypeSystem>::set_rtol(Scalar error) {
        rtol_ = error;
    }

    template <typename TypeSystem>
    template <typename Array1, typename Array2, typename Array3>
    auto RMS<TypeSystem>::error(const Array1 &y0, const Array2 &y1, const Array3 &y1_prime) -> Scalar {
        size_t const size = y0.size();
        Scalar error = 0;

        if constexpr (std::is_same_v<raw_type_t<typename Array1::value_type>, raw_type_t<Scalar>>) {
            for (size_t i = 0; i < size; ++i) {
                Scalar scale = std::max(fabs(y0[i]), fabs(y1[i])) * rtol_ + atol_;
                if (scale == 0) {
                    continue;
                }
                auto r = fabs(y1[i] - y1_prime[i]) / scale;
                error += r * r;
            }
        } else if constexpr (std::is_same_v<typename Array1::value_type, Vec3<Scalar>>) {
            for (size_t i = 0; i < size; ++i) {
                auto scale = vec_max(vec_abs(y0[i]), vec_abs(y1[i])) * rtol_ + atol_;
                if (scale == 0) {
                    continue;
                }
                auto v = vec_abs(y1[i] - y1_prime[i]) / scale;
                auto r = norm2(v) / 3;
                error += r;
                // Scalar scale = std::max(max_abs(y0[i]), max_abs(y1[i]));
                // auto r = max_abs(y1[i] - y1_prime[i]) / (atol_ + scale * rtol_);
                // error += r * r;
            }
        } else {
            spacehub_abort("Unsupported array type!");
        }
        return sqrt(error / size);
    }

    template <typename TypeSystem>
    RMS<TypeSystem>::RMS(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}
}  // namespace hub::ode
