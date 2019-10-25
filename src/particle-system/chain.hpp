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
 * @file chain.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_CHAIN_HPP
#define SPACEHUB_CHAIN_HPP

#include <algorithm>
#include <list>
#include <vector>
#include "../core-computation.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
      Class Chain Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 */
class Chain {
 public:
  struct Node {
    double r;   /**< Relative distance of two particles.*/
    size_t i;   /**< Particle index.*/
    size_t j;   /**< Particle index.*/
    bool avail; /**< State of node. If this node can be chained.*/
  };

  // Constructors
  SPACEHUB_MAKE_CONSTRUCTORS(Chain, default, default, default, default, default);

  template <typename Coord, typename IdxArray>
  static void calc_chain_index(Coord const &pos, IdxArray &index);

  template <typename Coord, typename IdxArray>
  static void update_chain(Coord &chain, IdxArray const &idx, IdxArray const &new_idx);

  template <typename ScalarArray, typename Coord, typename IdxArray>
  static void calc_cartesian(ScalarArray const &mass, Coord const &chain, Coord &cartesian, IdxArray const &index);

  template <typename Coord, typename IdxArray>
  static void calc_chain(Coord const &cartesian, Coord &chain, IdxArray const &index);

 private:
  static constexpr bool auto_centre_of_mass_move_{true};

  template <typename T>
  static bool not_in_list(std::list<T> &list, T var);

  template <typename InsertOpt>
  static bool try_insert(std::list<size_t> &list, size_t &chain_end, Node &n, size_t idx, InsertOpt insert);

  static bool try_add_to_chain(std::list<size_t> &list, size_t &head, size_t &tail, Node &n);

  template <typename Coord, typename Container>
  static void create_distances_array(Coord const &pos, Container &vec);

  template <typename IdxArray>
  static void create_index_from_dist_array(std::vector<Node> &dist, IdxArray &idx, size_t num);

  template <typename Coord>
  static auto get_new_node(Coord const &chain, size_t head, size_t tail) -> typename Coord::Vector;

  template <typename Array, typename IdxArray>
  static void to_chain(Array const &cartesian, Array &chain, IdxArray const &index);

  template <typename Array, typename IdxArray>
  static void to_cartesian(Array const &chain, Array &cartesian, IdxArray const &index);
};
/*---------------------------------------------------------------------------*\
      Class Chain Implementation
\*---------------------------------------------------------------------------*/
template <typename Coord, typename IdxArray>
void Chain::calc_chain_index(const Coord &pos, IdxArray &index) {
  std::vector<Node> dist;
  create_distances_array(pos, dist);
  std::sort(dist.begin(), dist.end(), [&](Node const &ni, Node const &nj) { return (ni.r < nj.r); });
  create_index_from_dist_array(dist, index, pos.size());
}

template <typename Coord, typename IdxArray>
void Chain::update_chain(Coord &chain, const IdxArray &idx, const IdxArray &new_idx) {
  using Vector = typename Coord::Vector;
  size_t size = chain.size();

  Coord new_chain;
  new_chain.reserve(size);

  auto get_idx = [&](auto var) -> auto { return std::find(idx.begin(), idx.end(), var) - idx.begin(); };

  for (size_t i = 0; i < size - 1; ++i) {
    auto first = get_idx(new_idx[i]);
    auto last = get_idx(new_idx[i + 1]);
    new_chain.emplace_back(get_new_node(chain, first, last));
  }

  if constexpr (!auto_centre_of_mass_move_) {
    new_chain.emplace_back(Vector(0, 0, 0));
  } else {
    Vector new_head = get_new_node(chain, 0, get_idx(new_idx[0]));
    new_chain.emplace_back(Vector(chain.x.back(), chain.y.back(), chain.z.back()) + new_head);
  }

  chain = std::move(new_chain);
}

template <typename ScalarArray, typename Coord, typename IdxArray>
void Chain::calc_cartesian(const ScalarArray &mass, const Coord &chain, Coord &cartesian, const IdxArray &index) {
  to_cartesian(chain.x, cartesian.x, index);
  to_cartesian(chain.y, cartesian.y, index);
  to_cartesian(chain.z, cartesian.z, index);
  if constexpr (!auto_centre_of_mass_move_) {
    calc::coord_move_to_com(mass, cartesian);
  }
}

template <typename Coord, typename IdxArray>
void Chain::calc_chain(const Coord &cartesian, Coord &chain, const IdxArray &index) {
  to_chain(cartesian.x, chain.x, index);
  to_chain(cartesian.y, chain.y, index);
  to_chain(cartesian.z, chain.z, index);
}

template <typename T>
bool Chain::not_in_list(std::list<T> &list, T var) {
  return list.end() == std::find(list.begin(), list.end(), var);
}

template <typename InsertOpt>
bool Chain::try_insert(std::list<size_t> &list, size_t &chain_end, Chain::Node &n, size_t idx, InsertOpt insert) {
  n.avail = false;
  if (not_in_list(list, idx)) {
    insert(idx);
    chain_end = idx;
    return true;
  } else {
    return false;
  }
}

bool Chain::try_add_to_chain(std::list<size_t> &list, size_t &head, size_t &tail, Chain::Node &n) {
  if (head == n.i) {
    return try_insert(list, head, n, n.j, [&](size_t idx) { list.emplace_front(idx); });
  } else if (head == n.j) {
    return try_insert(list, head, n, n.i, [&](size_t idx) { list.emplace_front(idx); });
  } else if (tail == n.i) {
    return try_insert(list, tail, n, n.j, [&](size_t idx) { list.emplace_back(idx); });
  } else if (tail == n.j) {
    return try_insert(list, tail, n, n.i, [&](size_t idx) { list.emplace_back(idx); });
  }
  return false;
}

template <typename Coord, typename Container>
void Chain::create_distances_array(const Coord &pos, Container &vec) {
  size_t num = pos.size();
  vec.reserve(num * (num - 1));
  for (size_t i = 0; i < num; ++i) {
    for (size_t j = i + 1; j < num; ++j) {
      auto dx = pos.x[j] - pos.x[i];
      auto dy = pos.y[j] - pos.y[i];
      auto dz = pos.z[j] - pos.z[i];
      vec.emplace_back(Node{dx * dx + dy * dy + dz * dz, i, j, true});
    }
  }
}

template <typename IdxArray>
void Chain::create_index_from_dist_array(std::vector<Node> &dist, IdxArray &idx, size_t num) {
  auto head = dist[0].i;
  auto tail = dist[0].j;
  dist[0].avail = false;

  std::list<size_t> idx_list{head, tail};

  size_t chained_num = 1;

  size_t dist_size = dist.size();
  for (size_t k = 1; k < dist_size; ++k) {
    if (dist[k].avail) {
      if (try_add_to_chain(idx_list, head, tail, dist[k])) {
        chained_num++;
        if (chained_num == num) {
          break;
        } else {
          k = 1;
        }
      }
    }
  }

  idx.clear();
  for (auto &l : idx_list) {
    idx.emplace_back(l);
  }
}

template <typename Coord>
auto Chain::get_new_node(const Coord &chain, size_t head, size_t tail) -> typename Coord::Vector {
  using Scalar = typename Coord::Scalar;
  using Vector = typename Coord::Vector;

  Scalar sign{1};

  if (head > tail) {
    std::swap(head, tail);
    sign = -1;
  }

  auto connect = [](auto const &array, auto first, auto last) -> auto {
    auto new_d = array[first];
    for (size_t j = first + 1; j < last; ++j) {
      new_d += array[j];
    }

    return new_d;
  };

  return Vector(sign * connect(chain.x, head, tail), sign * connect(chain.y, head, tail),
                sign * connect(chain.z, head, tail));
}

template <typename Array, typename IdxArray>
void Chain::to_chain(const Array &cartesian, Array &chain, const IdxArray &index) {
  const size_t size = cartesian.size();
  if constexpr (!auto_centre_of_mass_move_) {
    chain[size - 1] = 0;
  } else {
    chain[size - 1] = cartesian[index[0]];
  }
  for (size_t i = 0; i < size - 1; ++i) chain[i] = cartesian[index[i + 1]] - cartesian[index[i]];
}

template <typename Array, typename IdxArray>
void Chain::to_cartesian(const Array &chain, Array &cartesian, const IdxArray &index) {
  const size_t size = cartesian.size();
  if constexpr (!auto_centre_of_mass_move_) {
    cartesian[index[0]] = 0;
  } else {
    cartesian[index[0]] = chain[size - 1];
  }
  for (size_t i = 1; i < size; ++i) cartesian[index[i]] = cartesian[index[i - 1]] + chain[i - 1];
}
}  // namespace space
#endif
