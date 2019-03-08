#ifndef CHAIN_H
#define CHAIN_H

#include<vector>
#include<list>
#include<algorithm>

namespace SpaceH::chain {

/** @brief Struture to store the relative distance and index of two particles.*/
    struct Node {
        double r; /**< Relative distance of two particles.*/
        size_t i; /**< Particle index.*/
        size_t j; /**< Particle index.*/
        bool avail; /**< State of node. If this node can be chained.*/
    };

    template<typename Coord, typename Container>
    void create_distances_array(Coord const &pos, Container &vec) {
        size_t num = pos.size();
        vec.reserve(num * (num - 1));
        for (size_t i = 0; i < num; ++i) {
            for (size_t j = i + 1; j < num; ++j) {
                auto dx = pos.x[j] - pos.x[i];
                auto dy = pos.y[j] - pos.y[i];
                auto dz = pos.z[j] - pos.z[i];
                vec.emplace_back(dx * dx + dy * dy + dz * dz, i, j, true);
            }
        }
    }

    template<typename T>
    bool not_in_list(std::list<T> &list, T var) {
        return list.end() == std::find(list.begin(), list.end(), var)
    }

    template<typename InsertOpt>
    bool try_insert(std::list<size_t> &list, size_t &chain_end, Node &n, size_t idx, InsertOpt insert) {
        n.avail = false;
        if (not_in_list(list, idx)) {
            insert(idx);
            chain_end = idx;
            return true;
        } else {
            return false;
        }
    }

    bool try_add_to_chain(std::list<size_t> &list, size_t &head, size_t &tail, Node &n) {
        if (head == n.i) {
            return try_insert(list, head, n.j, [&](size_t idx) { list.emplace_front(idx); });
        } else if (head == n.j) {
            return try_insert(list, head, n.i, [&](size_t idx) { list.emplace_front(idx); });
        } else if (tail == n.i) {
            return try_insert(list, tail, n.j, [&](size_t idx) { list.emplace_back(idx); });
        } else if (tail == n.j) {
            return try_insert(list, tail, n.i, [&](size_t idx) { list.emplace_back(idx); });
        }
        return false;
    }

    templace<typename IdxArray>
    void create_index_from_dist_array(std::vector<Node> &dist, IdxArray &idx, size_t num) {
        auto head = dist[0].i;
        auto tail = dist[0].j;
        dist[0].avail = false;

        std::list<size_t> idx_list{head,tail};

        size_t chained_num = 1;

        size_t dist_size = dist.size();
        for (size_t k = 1; k < dist_size; ++k) {
            if (dist[k].avail) {
                if (try_add_to_chain(idx_list, head, tail, disk[k])) {
                    chained_num++;
                    if (chained_num == num) {
                        break;
                    } else {
                        k = 0;
                    }
                }
            }
        }

        idx.clear();
        for (auto &l : idx_list) {
            idx.emplace_back(l);
        }
    }

    template<typename Coord, typename IdxArray>
    void calc_chain_index(Coord const &pos, IdxArray &index) {
        std::vector<Node> dist;
        create_distances_array(pos, dist);
        std::sort(dist.begin(), dist.end(), [&](Node const &Ni, Node const &Nj) { return (Ni.r < Nj.r); });
        create_index_from_dist_array(dist, index, pos.size());
    }

    template<typename Coord, typename IdxArray, typename T>
    auto get_new_node(Coord &chain, IdxArray &idx, size_t head, size_t tail) -> typename Coord::Vector {
        using Scalar = typename Coord::Scalar;
        using Vector = typename Coord::Vector;

        Scalar sign{1};

        if (head > tail) {
            std::swap(head, tail);
            sign = -1;
        }

        auto connect = [](auto &array, auto first, auto last) -> auto {
            auto new_d = array[first];
            for (size_t j = first + 1; j < last; ++j)
                new_d += chain[j];
            return new_d;
        }

        return Vector(sign * connect(chain.x, head, tail), sign * connect(chain.y, head, tail), sign * connect(chain.z, head, tail))
    }

    template<typename Coord, typename IdxArray>
    void update_chain(Coord &chain, IdxArray &idx, IdxArray &new_idx) {
        using Vector = typename Coord::Vector;
        size_t size = chain.size();

        Coord new_chain;
        new_chain.reserve(size);

        auto get_idx = [&](auto var)->auto { return std::find(idx.begin(), idx.end(), var) - idx.begin(); };

        Vector new_head = get_new_node(chain, idx, 0, get_idx(new_idx[0]));

        new_chain.emplace_back(Vector(chain.x.back(), chain.y.back(), chain.z.back()) + new_head);

        for (size_t i = 0; i < size - 1; ++i) {
            auto first = get_idx(new_idx[i]);
            auto last  = get_idx(new_idx[i + 1]);
            new_chain.emplace_back(get_new_node(chain, idx, first, last));
        }

        chain = std::move(new_chain);
    }


    template<typename Array, typename IdxArray>
    void to_chain(Array const &cartesian, Array &chain, IdxArray const &index) {
        const size_t size = cartesian.size();
        //chain[N - 1] = 0;
        chain[size - 1] = cartesian[index[0]];
        for (size_t i = 0; i < size - 1; ++i)
            chain[i] = cartesian[index[i + 1]] - cartesian[index[i]];
    }

    template<typename Array, typename IdxArray>
    void to_cartesian(Array const &chain, Array &cartesian, IdxArray const &index) {
        const size_t size = cartesian.size();
        //data[index[0]] = 0;
        cartesian[index[0]] = chain[size - 1];
        for (size_t i = 1; i < size; ++i)
            cartesian[index[i]] = cartesian[index[i - 1]] + chain[i - 1];
    }

    template<typename Coord, typename IdxArray>
    void coord_calc_chain(Coord const &cartesian, Coord &chain, IdxArray const &index) {
        to_chain(cartesian.x, chain.x, index);
        to_chain(cartesian.y, chain.y, index);
        to_chain(cartesian.z, chain.z, index);
    }

    template<typename Coord, typename IdxArray>
    void coord_calc_cartesian(Coord const &chain, Coord &cartesian, IdxArray const &index) {
        to_cartesian(chain.x, cartesian.x, index);
        to_cartesian(chain.y, cartesian.y, index);
        to_cartesian(chain.z, cartesian.z, index);
    }


}
#endif

