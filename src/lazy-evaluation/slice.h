
#pragma once

#include <cstddef>

#include "lazy_expr.h"
namespace hub::lazy {

#define slice(BEGIN, END, ...) (Slice<BEGIN, END, ##__VA_ARGS__>())

#define EXPR_CREATE_SLICE_INTERFACE                                                                                  \
    template <int UNIQ(Begin), int UNIQ(End), int UNIQ(Stride) = 1>                                                  \
    inline constexpr auto operator[](const Slice<UNIQ(Begin), UNIQ(End), UNIQ(Stride)> &slice)                       \
        ->Slice_Expr<typename TYPE_OF_SELF, UNIQ(Begin), UNIQ(End), UNIQ(Stride), TYPE_OF_SELF::size()> {            \
        static_assert(UNIQ(Begin) <= static_cast<int>(TYPE_OF_SELF::size()), "Slice out of range");                  \
        static_assert(UNIQ(Begin) >= -static_cast<int>(TYPE_OF_SELF::size()), "Slice out of range");                 \
        static_assert(UNIQ(End) <= static_cast<int>(TYPE_OF_SELF::size()), "Slice out of range");                    \
        static_assert(UNIQ(End) >= -static_cast<int>(TYPE_OF_SELF::size()), "Slice out of range");                   \
        return Slice_Expr<typename TYPE_OF_SELF, UNIQ(Begin), UNIQ(End), UNIQ(Stride), TYPE_OF_SELF::size()>(*this); \
    }

    constexpr size_t constexpr_slice_len(int begin, int end, int stride, size_t size) {
        const int b = begin > 0 ? begin : size + begin;
        const int e = end > 0 ? end : size + end;
        // return static_cast<size_t>(( (b < e) ? (e - b) : (b - e) ) / ((stride > 0) ? stride : -stride));*/
        if (stride > 0 && e > b) {
            return 1 + static_cast<size_t>(e - b);
        } else if (stride > 0 && e < b) {
            return 1 + static_cast<size_t>(e - b + size);
        } else if (stride < 0 && e > b) {
            return 1 + static_cast<size_t>(size + b - e);
        } else if (stride < 0 && e < b) {
            return 1 + static_cast<size_t>(b - e);
        } else {
            return 1;
        }
    }

    template <int Begin, int End, int Stride = 1>
    struct Slice {};

    template <typename T, int Begin, int End, int Stride, size_t SrcSize>
    struct Slice_Expr : public Expr<Slice_Expr<T, Begin, End, Stride, SrcSize>> {
       private:
        T &src_;

       public:
        using value_type = typename T::value_type;

        explicit Slice_Expr(T &src) : src_(src) {}

        constexpr static size_t size() { return constexpr_slice_len(Begin, End, Stride, SrcSize); }

        Slice_Expr &operator=(const value_type &scalar) {
            for (size_t i = 0; i < size(); ++i) {
                (*this)[i] = scalar;
            }
            return *this;
        }

        inline value_type const &eval(size_t i) const { return src_.eval((Begin + SrcSize + i * Stride) % SrcSize); }

        inline value_type const &operator[](size_t i) const { return eval(i); };

        inline value_type &operator[](size_t i) { return const_cast<value_type &>(eval(i)); };

        T max() const {
            T max_value = (*this)[0];
            for (size_t i = 1; i < size(); ++i) {
                if ((*this)[i] > max_value) max_value = (*this)[i];
            }
            return max_value;
        }

        T min() const {
            T min_value = (*this)[0];
            for (size_t i = 1; i < size(); ++i) {
                if ((*this)[i] < min_value) min_value = (*this)[i];
            }
            return min_value;
        }

        EXPR_CREATE_SLICE_INTERFACE;  // create operator[slice(begin,end,stride)]

        template <typename U>
        Slice_Expr &operator=(const Expr<U> &rhs_expr) {
            const U &rhs = rhs_expr.cast();
            multi_threads_loop(size(), multi_thread::machine_thread_num, [&](size_t begin, size_t end) {
                for (size_t i = begin; i < end; ++i) (*this)[i] = rhs.eval(i);
            });
            return *this;
        }
    };
}  // namespace hub::lazy
