#ifndef TYPECLASS_H
#define TYPECLASS_H

#include "vector/vector3.h"
#include "dev-tools.h"
#include <vector>

namespace space {

    template<typename T>
    struct Coords {
        using Scalar = typename T::value_type;
        using Vector = Vec3<Scalar>;

        Coords() = default;

        explicit Coords(size_t count) : x(count), y(count), z(count) {}

        Coords(Coords const &other) = default;

        Coords(Coords &&other) noexcept = default;

        Coords &operator=(Coords const &other) = default;

        Coords &operator=(Coords &&other) noexcept = default;

        size_t size() const {
            return x.size();
        }

        void reserve(size_t new_cap) {
            reserve_all(new_cap, x, y, z);
        }

        void resize(size_t new_sz) {
            resize_all(new_sz, x, y, z);
        }

        template<typename Vector>
        void emplace_back(Vector const &v) {
            x.emplace_back(v.x);
            y.emplace_back(v.y);
            z.emplace_back(v.z);
        }

        void emplace_back(Scalar &&xx, Scalar &&yy, Scalar &&zz) {
            x.emplace_back(std::forward<Scalar>(xx));
            y.emplace_back(std::forward<Scalar>(yy));
            z.emplace_back(std::forward<Scalar>(zz));
        }

        void shrink_to_fit() {
            shrink_to_fit_all(x, y, z);
        }

        void clear() {
            clear_all(x, y, z);
        }

        T x;
        T y;
        T z;
    };

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
    inline auto distance(Coords<T> const &c, size_t i, size_t j) {
        auto dx = c.x[i] - c.x[j];
        auto dy = c.y[i] - c.y[j];
        auto dz = c.z[i] - c.z[j];
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    /**
     * @brief Type system that is used across the whole spacehub system. This type class provide all basic type i.e
     * 'Scalar', 'Vector', 'ScalarArray', 'Coord' and etc,.
     *
     * @tparam Real The basic scalar type.
     * @tparam TContainer The container type.
     */
    template<typename Real, template<class...> class TContainer = std::vector>
    struct Types {
    public:
        template<typename ...T>
        using Container = TContainer<T...>;

        using Scalar      = Real;
        using ScalarArray = Container<Scalar>;
        using IntArray    = Container<int>;
        using IdxArray    = Container<size_t>;
        using Vector      = Vec3<Scalar>;
        using VectorArray = Container<Vector>;
        using Coord       = Coords<ScalarArray>;
    };
}//end namespace SpaceH

#endif
