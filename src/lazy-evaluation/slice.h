
#ifndef SPACEHUB_SLICE_H
#define SPACEHUB_SLICE_H

#include <cstddef>
#include "lazy_expr.h"

namespace SpaceH {
    namespace Lazy {

#define slice(BEGIN, END, ...) (Slice<BEGIN,END ,##__VA_ARGS__>())

#define EXPR_CREATE_SLICE_INTERFACE                                                                                    \
    template<int UNIQ(Begin), int UNIQ(End), int UNIQ(Stride) = 1>                                                     \
    inline constexpr auto operator[](const Slice<UNIQ(Begin), UNIQ(End), UNIQ(Stride)>& slice)                         \
    ->Slice_Expr<typename TYPE_OF_SELF, UNIQ(Begin), UNIQ(End), UNIQ(Stride), TYPE_OF_SELF::len>{                      \
        static_assert(UNIQ(Begin) <= static_cast<int>(TYPE_OF_SELF::len), "Slice out of range");                       \
        static_assert(UNIQ(Begin) >= -static_cast<int>(TYPE_OF_SELF::len), "Slice out of range");                      \
        static_assert(UNIQ(End)   <= static_cast<int>(TYPE_OF_SELF::len), "Slice out of range");                       \
        static_assert(UNIQ(End)   >= -static_cast<int>(TYPE_OF_SELF::len), "Slice out of range");                      \
        return Slice_Expr<typename TYPE_OF_SELF, UNIQ(Begin), UNIQ(End), UNIQ(Stride), TYPE_OF_SELF::len>(*this);      \
    }

        constexpr size_t constexpr_slice_len(int begin, int end, int stride, size_t size){
            const int b = begin > 0 ? begin : size + begin;
            const int e = end > 0 ? end : size + end;
            //return static_cast<size_t>(( (b < e) ? (e - b) : (b - e) ) / ((stride > 0) ? stride : -stride));*/
            if (stride > 0 && e > b){
                return 1 + static_cast<size_t>(e - b);
            } else if(stride > 0 && e < b){
                return 1 + static_cast<size_t>(e - b + size);
            } else if(stride < 0 && e > b){
                return 1 + static_cast<size_t>(size + b - e);
            } else if(stride < 0 && e < b){
                return 1 + static_cast<size_t>(b - e);
            } else {
                return 1;
            }
        }

        template<int Begin, int End, int Stride = 1>
        struct Slice{};

        template <typename T, int Begin, int End, int Stride, size_t SrcSize>
        struct Slice_Expr : public Expr<Slice_Expr<T, Begin, End, Stride, SrcSize>>{
        private:
            T& src_;
        public:
            static constexpr size_t len{constexpr_slice_len(Begin,End, Stride, SrcSize)};
            using value_type = typename  T::value_type;

            Slice_Expr(T& src)
                : src_(src) {}

            const Slice_Expr& operator=(const value_type& single){
                for (size_t i = 0; i < len; ++i) {
                    (*this).dst(i) = single;
                }
                return *this;
            }
            inline const value_type& eval(size_t i) const {
                return src_.eval((Begin + SrcSize + i*Stride)%SrcSize);
            }
            inline value_type& dst(size_t i) {
                return src_.dst((Begin + SrcSize + i*Stride)%SrcSize);
            };

            EXPR_CREATE_SLICE_INTERFACE;//create operator[slice(begin,end,stride)]
            EXPR_CREATE_LINEAR_ASSIGN_OPERATORS;//create operator=,+=,-=,*=,/= for Expr<T> as rhs.
        };

    }
}
#endif //SPACEHUB_SLICE_H
