//
// Created by yihan on 4/25/19.
//

#ifndef SPACEHUB_COORDS_H
#define SPACEHUB_COORDS_H

#include "vector/vector3.h"

namespace space {

    /*---------------------------------------------------------------------------*\
        Class Coords Declaration
    \*---------------------------------------------------------------------------*/
    template<typename T>
    struct Coords {
        //type members
        using Scalar = typename T::value_type;
        using Vector = Vec3<Scalar>;

        //constructors
        Coords() = default;

        explicit Coords(size_t count);

        Coords(Coords const &other) = default;

        Coords(Coords &&other) noexcept = default;

        Coords &operator=(Coords const &other) = default;

        Coords &operator=(Coords &&other) noexcept = default;


        //public methods
        void clear();

        template<typename GenVector>
        void emplace_back(GenVector const &v);

        void emplace_back(Scalar &&xx, Scalar &&yy, Scalar &&zz);

        void reserve(size_t new_cap);

        void resize(size_t new_sz);

        void shrink_to_fit();

        size_t size() const;

        //public members
        T x;
        T y;
        T z;
    };

    /*---------------------------------------------------------------------------*\
        Class Coords Implementation
    \*---------------------------------------------------------------------------*/

    template<typename T>
    Coords<T>::Coords(size_t count) : x(count), y(count), z(count) {}

    template<typename T>
    size_t Coords<T>::size() const {
        return x.size();
    }

    template<typename T>
    void Coords<T>::reserve(size_t new_cap) {
        reserve_all(new_cap, x, y, z);
    }

    template<typename T>
    void Coords<T>::resize(size_t new_sz) {
        resize_all(new_sz, x, y, z);
    }

    template<typename T>
    template<typename GenVector>
    void Coords<T>::emplace_back(GenVector const &v) {
        x.emplace_back(v.x);
        y.emplace_back(v.y);
        z.emplace_back(v.z);
    }

    template<typename T>
    void Coords<T>::emplace_back(Scalar &&xx, Scalar &&yy, Scalar &&zz) {
        x.emplace_back(std::forward<Scalar>(xx));
        y.emplace_back(std::forward<Scalar>(yy));
        z.emplace_back(std::forward<Scalar>(zz));
    }

    template<typename T>
    void Coords<T>::shrink_to_fit() {
        shrink_to_fit_all(x, y, z);
    }

    template<typename T>
    void Coords<T>::clear() {
        clear_all(x, y, z);
    }

    /*---------------------------------------------------------------------------*\
        Help functions and tools
    \*---------------------------------------------------------------------------*/
    template<typename STL, typename T>
    void add_coords_to(STL &stl, Coords<T> const &coords) {
        for (auto const &xx : coords.x) {
            stl.emplace_back(xx);
        }
        for (auto const &yy : coords.y) {
            stl.emplace_back(yy);
        }
        for (auto const &zz : coords.z) {
            stl.emplace_back(zz);
        }
    }

    template<typename STLIterator, typename T>
    void load_to_coords(STLIterator &i, Coords<T> &coords) {
        for (auto &xx : coords.x) {
            xx = *i;
            i++;
        }
        for (auto &yy : coords.y) {
            yy = *i;
            i++;
        }
        for (auto &zz : coords.z) {
            zz = *i;
            i++;
        }
    }

    template<typename T>
    auto distance(Coords<T> const &c, size_t i, size_t j) {
        auto dx = c.x[i] - c.x[j];
        auto dy = c.y[i] - c.y[j];
        auto dz = c.z[i] - c.z[j];
        return sqrt(dx * dx + dy * dy + dz * dz);
    }
}

#endif //SPACEHUB_COORDS_H
