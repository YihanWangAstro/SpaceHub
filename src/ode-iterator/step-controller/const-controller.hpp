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
 * @file const-controller.hpp
 *
 * Header file.
 */
#pragma once

#include "../../math.hpp"

namespace space::ode_iterator {
    /*---------------------------------------------------------------------------*\
         Class Const step Controller Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief
     *
     * @tparam TypeSystem
     */
    template <typename TypeSystem>
    class ConstStepController {
       public:
        // Type member
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        template <typename ArrayLike>
        inline Scalar next_with_limiter(size_t order, Scalar old_step, ArrayLike const &errors) {
            return old_step;
        };

        inline Scalar next_with_limiter(size_t order, Scalar old_step, Scalar error) { return old_step; };

        template <typename ArrayLike>
        inline Scalar next(size_t order, Scalar old_step, ArrayLike const &errors) {
            return old_step;
        };

        inline Scalar next(size_t order, Scalar old_step, Scalar error) { return old_step; };
    };
}  // namespace space::ode_iterator
