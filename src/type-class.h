#ifndef TYPECLASS_H
#define TYPECLASS_H

#include "vector/vector3.h"
#include "dev-tools.h"
#include <array>
#include <vector>

namespace SpaceH {

    template<typename T>
    struct Coords{
        using Scalar = typename T::value_type;
        using Vector = Vec3<Scalar>;

        Coords() = default;

        explicit Coords(size_t count) : x(count), y(count), z(count) {}

        Coords(Coords const & other) : x(other.x), y(other.y), z(other.z) {}

        Coords(Coords&& other) noexcept: x(std::move(other.x)), y(std::move(other.x)), z(std::move(other.z)) {}

        Coords& operator=(Coords const& other) {
            x = other.x;
            x = other.y;
            x = other.z;
            return *this;
        }

        Coords& operator=(Coords && other) noexcept {
            x = std::move(other.x);
            x = std::move(other.y);
            x = std::move(other.z);
            return *this;
        }

        size_t size() const {
            return x.size();
        }

        void reserve(size_t new_cap){
            x.reserve(new_cap);
            y.reserve(new_cap);
            z.reserve(new_cap);
        }

        void resize(size_t new_sz){
            x.resize(new_sz);
            y.resize(new_sz);
            z.resize(new_sz);
        }

        template<typename Vector>
        void emplace_back(Vector const& v){
            x.emplace_back(v.x);
            y.emplace_back(v.y);
            z.emplace_back(v.z);
        }

        void emplace_back(Scalar && xx, Scalar && yy, Scalar && zz){
            x.emplace_back(std::forward<Scalar>(xx));
            y.emplace_back(std::forward<Scalar>(yy));
            z.emplace_back(std::forward<Scalar>(zz));
        }

        void shrink_to_fit() {
            x.shrink_to_fit();
            y.shrink_to_fit();
            z.shrink_to_fit();
        }

        void clear() {
            x.clear();
            y.clear();
            z.clear();
        }

        T x;
        T y;
        T z;
    };

    template <typename T>
    inline auto distance(Coords<T> const& c, size_t i, size_t j){
        auto dx = c.x[i] - c.x[j];
        auto dy = c.y[i] - c.y[j];
        auto dz = c.z[i] - c.z[j];
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

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
