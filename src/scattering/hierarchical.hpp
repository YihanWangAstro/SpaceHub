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
 * @file hierarchical.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_HIERARCHICAL_HPP
#define SPACEHUB_HIERARCHICAL_HPP

#include "../orbits/orbits.hpp"
#include "../vector/vector3.hpp"

/**
 * @namespace space::scattering
 * namespace for scattering
 */
namespace space::scattering {

template <typename Vector>
class HierarchicalNode {
 public:
  Vector pos;
  Vector vel;
  double mass;
  std::string name;
};
/*template <typename Node>
void node_swap(Node& n1, Node& n2) {
  std::swap(n1.mass, n2.mass);
  std::swap(n1.name, n2.name);
}*/
template <typename Particles, typename Node>
void create_init_list(Particles const& ptc, std::vector<Node>& vec) {
  size_t p_num = ptc.number();
  vec.clear();
  vec.reserve(p_num);
  for (size_t i = 0; i < p_num; ++i) {
    vec.emplace_back(Node{ptc[i].pos, ptc[i].vel, ptc[i].m, std::to_string(i)});
  }
}

template <typename Nodes>
bool check_most_bound(Nodes const& vec, size_t idx, double amin) {
  for (size_t i = 1; i < vec.size(); i++) {
    if (i != idx) {
      auto [a, e] = orbit::calc_a_e(vec[idx].mass + vec[i].mass, vec[idx].pos - vec[i].pos, vec[idx].vel - vec[i].vel);
      if (0 < a && a < amin) {
        return false;
      }
    }
  }
  return true;
}

template <typename Particles>
std::string to_hierarchical(Particles const& ptc) {
  using Vector = decltype(ptc[0].pos);
  using Node = HierarchicalNode<Vector>;
  std::vector<Node> vec{0};
  create_init_list(ptc, vec);
  std::string out;
  for (; vec.size() > 0;) {
    Vector p0 = vec[0].pos;
    Vector v0 = vec[0].vel;
    double m0 = vec[0].mass;
    double amin = 1e999;
    size_t idx = 0;
    for (size_t i = 1; i < vec.size(); i++) {
      auto [a, e] = orbit::calc_a_e(m0 + vec[i].mass, p0 - vec[i].pos, v0 - vec[i].vel);
      if (0 < a && a < amin) {
        idx = i;
        amin = a;
      }
    }

    if (idx == 0) {
      out += vec[0].name;
      vec.erase(vec.begin());
    } else {
      if (check_most_bound(vec, idx, amin)) {
        double m1 = vec[idx].mass;
        double mt = m1 + m0;
        Vector p1 = vec[idx].pos;
        Vector v1 = vec[idx].vel;
        vec.emplace_back(Node{(m0 * p0 + m1 * p1) / mt, (m0 * v0 + m1 * p1) / mt, mt,
                              "(" + vec[0].name + "," + vec[idx].name + ")"});
        vec.erase(vec.begin() + idx);
        vec.erase(vec.begin());
      } else {
        std::swap(vec[0], vec[idx]);
      }
    }
  }
  return out;
}

}  // namespace space::scattering

#endif  // SPACEHUB_CROSS_SECTION_HPP
