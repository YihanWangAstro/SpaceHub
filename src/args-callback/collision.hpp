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
 *
 * @file collision.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_COLLISION_HPP
#define SPACEHUB_COLLISION_HPP
#include "../dev-tools.hpp"
#include "../vector/vector3.hpp"
namespace run_operations {
CREATE_METHOD_CHECK(radius);
class StickyCollision {
  template <typename ParticleSys>
  auto operator()(ParticleSys const& ptc) {
    static_assert(HAS_METHOD(PartileSys, radius), "Collision detection needs the particles have the radius!");
    auto const& p = ptc.pos();
    auto const& r = ptc.radius();
    size_t num = ptc.number();
    for (size_t i = 0 ; i < num; ++i){
      for (size_t j  = i+1; j < num ; ++j){
        if( distance(p[j],p[i]) < r[i] + r[j]) 
          return true;
      }
    } 
return false;
    }
};
}  // namespace run_operations
#endif  // SPACEHUB_COLLISION_HPP
