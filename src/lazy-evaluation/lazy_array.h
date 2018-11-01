
#ifndef SPACEHUB_LAZY_ARRAY_H
#define SPACEHUB_LAZY_ARRAY_H

#include "lazy_expr.h"
#include "slice.h"

namespace SpaceH {
    namespace Lazy {

        template<typename T>
        constexpr bool leq_cache_line(size_t len = 1) {
            return sizeof(T) * len <= sizeof(double) * 4;
        }

        template<typename Element, size_t Len, bool IsSmall = leq_cache_line<Element>(Len)>
        struct Larray : public Expr<Larray<Element, Len, IsSmall>> {
        public:
            static constexpr size_t len{Len};
            using value_type = Element;

            Larray()
                    : data_(new Element[Len]) {}

            Larray(const Larray &src)
                    : data_(new Element[Len]) {
                for (size_t i = 0; i < len; ++i) {
                    data_[i] = src[i];
                }
            }

            Larray(Larray &&src)
                    : data_(std::move(src.data_)) {}

            template<typename T>
            Larray(const Expr<T> &expr)
                    : data_(new Element[Len]) {
                const T &src = expr.cast();
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = src[i];
                }
            }

            template<size_t S>
            Larray(const Element (&src)[S])
                    : data_(new Element[Len]) {
                COMPILE_TIME_ASSERT(S == Len, "Size inconsistency of initializer and array!");
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = src[i];
                }
            }

            Element const & eval(size_t i) const {
                return data_[i];
            }

            Element& dst(size_t i){
                return data_[i];
            }

            template<typename... Args>
            Larray(Element arg, Args ... args)
                    : data_(new Element[Len]{arg, static_cast<Element>(args) ...}) {}

            inline static constexpr size_t size() {
                return Len;
            }

            inline const Element &operator[](size_t i) const {
                return data_[i];
            }

            inline Element &operator[](size_t i) {
                return const_cast<Element&>(static_cast<Larray const&>(*this)[i]);
            }

            const Larray &operator=(const Larray &src) {
                for (size_t i = 0; i < len; ++i) {
                    data_[i] = src[i];
                }
                return *this;
            }

            const Larray &operator=(Larray &&src) {
                data_.reset();
                data_ = std::move(src.data_);
                return *this;
            }

            EXPR_CREATE_SLICE_INTERFACE;//create operator[slice(begin,end,stride)]
            EXPR_CREATE_LINEAR_ASSIGN_OPERATORS;//create operator=,+=,-=,*=,/= for Expr<T> as rhs.
        private:
            std::unique_ptr<Element[]> data_;
        };


        template<typename Element, size_t Len>
        struct Larray<Element, Len, true> : public Expr<Larray<Element, Len, true>> {
        public:
            static constexpr size_t len{Len};
            using value_type = Element;

            Larray() = default;

            Larray(const Larray &src) = default;

            Larray(const Element &single) {
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = single;
                }
            }

            template<size_t S>
            Larray(const Element (&src)[S]) {
                COMPILE_TIME_ASSERT(S == Len, "Size inconsistency of initializer and array!");
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = src[i];
                }
            }

            template<typename T>
            Larray(const Expr<T> &expr) {
                const T &src = expr.cast();
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = src[i];
                }
            }

            template<typename... Args>
            Larray(Element arg, Args ... args)
                    : data_{arg, static_cast<Element>(args) ...} {}

            Element const& eval(size_t i){
                return data_[i];
            }

            Element& dst(size_t i){
                return data_[i];
            }

            inline static constexpr size_t size() {
                return Len;
            }

            inline Element const &operator[](size_t i) const {
                return data_[i];
            }

            inline Element &operator[](size_t i) {
                return const_cast<Element&>(static_cast<Larray const&>(*this)[i]);
            }
            EXPR_CREATE_SLICE_INTERFACE;//create operator[slice(begin,end,stride)]
            EXPR_CREATE_LINEAR_ASSIGN_OPERATORS;//create operator=,+=,-=,*=,/= for Expr<T> as rhs.
        private:
            Element data_[Len];
        };
    }
}

#endif //SPACEHUB_LAZY_ARRAY_H
