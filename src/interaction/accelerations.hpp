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
 * @file accelerations.hpp
 * Header file.
 */
#include "../dev-tools.hpp"

#ifndef SPACEHUB_ACCELERATIONS_H
#define SPACEHUB_ACCELERATIONS_H
namespace space::interactions {
/*---------------------------------------------------------------------------*\
    Class Accelerations Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Interactions
 * @tparam Coord
 */
template <typename Interactions, typename Coord>
class Accelerations {
 public:
  // Constructors
  SPACEHUB_MAKE_CONSTRUCTORS(Accelerations, default, default, default, default, default);

  explicit Accelerations(size_t size);

  // Public methods
  SPACEHUB_STD_ACCESSOR(auto, acc, acc_);

  SPACEHUB_STD_ACCESSOR(auto, newtonian_acc, newtonian_acc_);

  SPACEHUB_STD_ACCESSOR(auto, tot_vel_indep_acc, tot_vel_indep_acc_);

  SPACEHUB_STD_ACCESSOR(auto, ext_vel_indep_acc, ext_vel_indep_acc_);

  SPACEHUB_STD_ACCESSOR(auto, ext_vel_dep_acc, ext_vel_dep_acc_);

 private:
  Coord acc_;

  Coord newtonian_acc_;

  Coord tot_vel_indep_acc_;

  std::conditional_t<Interactions::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;

  std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;
};

template <typename Interactions, typename Coord>
Accelerations<Interactions, Coord>::Accelerations(size_t size)
    : acc_{size}, newtonian_acc_{size}, tot_vel_indep_acc_{size} {
  if constexpr (Interactions::ext_vel_indep) {
    ext_vel_indep_acc_.resize(size);
  }
  if constexpr (Interactions::ext_vel_dep) {
    ext_vel_dep_acc_.resize(size);
  }
}

/*---------------------------------------------------------------------------*\
    Class Accelerations Implementation
\*---------------------------------------------------------------------------*/
}  // namespace space::interactions
#endif  // SPACEHUB_ACCELERATIONS_H
