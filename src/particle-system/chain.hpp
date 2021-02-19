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
 * @file chain.hpp
 *
 * Header file.
 */
#pragma once

#include <algorithm>
#include <list>
#include <vector>

#include "../core-computation.hpp"

namespace space {

    /*---------------------------------------------------------------------------*\
          Class Chain Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief Static class for chain coordinates transformation.
     *
     */
    class Chain {
       public:
        /**
         * @brief Node used to create chain index.
         */
        struct Node {
            /**
             * @brief Relative distance between two particles.
             */
            double r;
            /**
             * @brief The index of the first particle.
             */
            size_t i;
            /**
             * @brief The index of the second particle.
             */
            size_t j;
            /**
             * @brief Flag indicates if this node has been added to the chain.
             */
            bool avail;
        };

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(Chain, default, default, default, default, default);

        /**
         * @brief Calculate the chain index.
         *
         * @tparam VectorArray Type of the Structure of Array coordinates.
         * @tparam IdxArray Type of the index array.
         * @param[in] pos Input position in Cartesian coordinates.
         * @param[out] index Output index array.
         * @note The memory space of the index need to be allocated in advance. The size of pos and index need to be the
         * same.
         */
        template <typename VectorArray, typename IdxArray>
        static void calc_chain_index(VectorArray const &pos, IdxArray &index);

        /**
         * @brief Update the chain coordinates from old index array to new index array.
         *
         * @tparam VectorArray Type of the Structure of Array coordinates.
         * @tparam IdxArray Type of the index array.
         * @param[in,out] chain The chain coordinates.
         * @param[in] idx Old chain index array.
         * @param[in] new_idx New chain index array.
         */
        template <typename VectorArray, typename IdxArray>
        static void update_chain(VectorArray &chain, VectorArray const &cartesian, IdxArray const &idx,
                                 IdxArray const &new_idx);

        /**
         * @brief Calculate the corresponding Cartesian coordinates of the chain coordinates.
         *
         * @tparam ScalarArray Type of the mass array.
         * @tparam VectorArray Type of the Structure of Array coordinates.
         * @tparam IdxArray Type of the index array.
         * @param[in] mass The mass array(needed for centre of mass movement).
         * @param[in] chain The chain coordinates.
         * @param[out] cartesian The Cartesian coordinates.
         * @param[in] index The chain index array.
         * @note The memory space of the cartesian need to be allocated in advance. The size of the input parameters
         * need to be the same.
         */
        template <typename ScalarArray, typename VectorArray, typename IdxArray>
        static void calc_cartesian(ScalarArray const &mass, VectorArray const &chain, VectorArray &cartesian,
                                   IdxArray const &index);

        /**
         * @brief Calculate the corresponding chain coordinates of the Cartesian coordinates.
         *
         * @tparam VectorArray Type of the Structure of Array coordinates.
         * @tparam IdxArray Type of the index array.
         * @param[in] cartesian The Cartesian coordinates.
         * @param[out] chain The chain coordinates.
         * @param[in] index The chain index array.
         */
        template <typename VectorArray, typename IdxArray>
        static void calc_chain(VectorArray const &cartesian, VectorArray &chain, IdxArray const &index);

        /**
         * @brief Method flag. Indicates the transformation type;
         *
         * If this flag is true, the transformation between chain coordinates and Cartesian coordinates is bijective
         * mapping, otherwise, a centre of mass movement will be performed after the transformation from chain
         * coordinates to Cartesian coordinates.
         */
        static constexpr bool bijective_transfer{true};

       private:
        template <typename T>
        static bool not_in_list(std::list<T> &list, T var);

        template <typename InsertOpt>
        static bool try_insert(std::list<size_t> &list, size_t &chain_end, Node &n, size_t idx, InsertOpt insert);

        static bool try_add_to_chain(std::list<size_t> &list, size_t &head, size_t &tail, Node &n);

        template <typename VectorArray, typename Container>
        static void create_distances_array(VectorArray const &pos, Container &vec);

        template <typename IdxArray>
        static void create_index_from_dist_array(std::vector<Node> &dist, IdxArray &idx, size_t num);

        template <typename VectorArray>
        static auto get_new_node(VectorArray const &chain, size_t head, size_t tail) ->
            typename VectorArray::value_type;

        template <typename Array, typename IdxArray>
        static void to_chain(Array const &cartesian, Array &chain, IdxArray const &index);

        template <typename Array, typename IdxArray>
        static void to_cartesian(Array const &chain, Array &cartesian, IdxArray const &index);

        CREATE_MEMBER_CHECK(err);
    };

    /*---------------------------------------------------------------------------*\
          Class Chain Implementation
    \*---------------------------------------------------------------------------*/
    template <typename VectorArray, typename IdxArray>
    void Chain::calc_chain_index(VectorArray const &pos, IdxArray &index) {
        std::vector<Node> dist;
        create_distances_array(pos, dist);
        std::sort(dist.begin(), dist.end(), [&](Node const &ni, Node const &nj) { return (ni.r < nj.r); });
        create_index_from_dist_array(dist, index, pos.size());
    }

    template <typename VectorArray, typename IdxArray>
    void Chain::update_chain(VectorArray &chain, VectorArray const &cartesian, const IdxArray &idx,
                             const IdxArray &new_idx) {
        using Vector = typename VectorArray::value_type;
        using Scalar = typename Vector::value_type;

        size_t size = chain.size();

        VectorArray new_chain;
        new_chain.reserve(size);

        auto get_idx = [&](auto var) -> auto { return std::find(idx.begin(), idx.end(), var) - idx.begin(); };

        for (size_t i = 0; i < size - 1; ++i) {
            auto first = get_idx(new_idx[i]);
            auto last = get_idx(new_idx[i + 1]);
            new_chain.emplace_back(get_new_node(chain, first, last));
        }

        if constexpr (!bijective_transfer) {
            new_chain.emplace_back(Vector(0, 0, 0));
        } else {
            new_chain.emplace_back(cartesian[new_idx[0]]);
        }
        chain = std::move(new_chain);
    }

    template <typename ScalarArray, typename VectorArray, typename IdxArray>
    void Chain::calc_cartesian(const ScalarArray &mass, VectorArray const &chain, VectorArray &cartesian,
                               const IdxArray &index) {
        to_cartesian(chain, cartesian, index);
        if constexpr (!bijective_transfer) {
            calc::move_to_com(mass, cartesian);
        }
    }

    template <typename VectorArray, typename IdxArray>
    void Chain::calc_chain(VectorArray const &cartesian, VectorArray &chain, const IdxArray &index) {
        to_chain(cartesian, chain, index);
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

    template <typename VectorArray, typename Container>
    void Chain::create_distances_array(VectorArray const &pos, Container &vec) {
        size_t num = pos.size();
        vec.reserve(num * (num - 1));
        for (size_t i = 0; i < num; ++i) {
            for (size_t j = i + 1; j < num; ++j) {
                decltype(pos[j]) dr = (pos[j] - pos[i]);
                decltype(dr.x) r2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                vec.emplace_back(Node{r2, i, j, true});
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

    template <typename VectorArray>
    auto Chain::get_new_node(VectorArray const &chain, size_t head, size_t tail) -> typename VectorArray::value_type {
        using Vector = typename VectorArray::value_type;
        using Scalar = typename Vector::value_type;

        bool sign{true};

        if (head > tail) {
            std::swap(head, tail);
            sign = false;
        }

        auto connect = [](auto const &array, auto first, auto last, bool sign) -> auto {
            if (sign) {
                Vector new_d = array[first];
                if constexpr (HAS_MEMBER(typename Vector::value_type, err)) {
                    new_d.x.err = new_d.y.err = new_d.z.err = 0;
                }
                for (size_t j = first + 1; j < last; ++j) {
                    new_d = new_d + array[j];
                    // new_d += array[j];
                }
                return new_d;
            } else {
                Vector new_d = -array[first];
                if constexpr (HAS_MEMBER(typename Vector::value_type, err)) {
                    new_d.x.err = new_d.y.err = new_d.z.err = 0;
                }

                for (size_t j = first + 1; j < last; ++j) {
                    new_d = new_d - array[j];
                    // new_d -= array[j];
                }
                return new_d;
            }
        };

        return connect(chain, head, tail, sign);
    }

    template <typename Array, typename IdxArray>
    void Chain::to_chain(const Array &cartesian, Array &chain, const IdxArray &index) {
        using T = typename Array::value_type;
        const size_t size = cartesian.size();
        if constexpr (!bijective_transfer) {
            chain[size - 1] = T{0.0};
        } else {
            chain[size - 1] = cartesian[index[0]];
        }
#pragma GCC ivdep
        for (size_t i = 0; i < size - 1; ++i) {
            chain[i] = cartesian[index[i + 1]] - cartesian[index[i]];
        }
    }

    template <typename Array, typename IdxArray>
    void Chain::to_cartesian(const Array &chain, Array &cartesian, const IdxArray &index) {
        using T = typename Array::value_type;
        const size_t size = cartesian.size();
        if constexpr (!bijective_transfer) {
            cartesian[index[0]] = T{0.0};
        } else {
            cartesian[index[0]] = chain[size - 1];
        }
        for (size_t i = 1; i < size; ++i) {
            cartesian[index[i]] = cartesian[index[i - 1]] + chain[i - 1];
        }
    }
}  // namespace space
