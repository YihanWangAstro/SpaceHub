
#ifndef SPACEHUB_LAZY_EXPRESSION_H
#define SPACEHUB_LAZY_EXPRESSION_H

#include "dev_tools.h"
#include "stddef.h"
#include <math.h>

namespace SpaceH {
    namespace Lazy {

        template<typename Lhs, typename Rhs>
        constexpr size_t constexpr_size(){
            if constexpr (INDEXABLE(Lhs)){
                return Lhs::size;
            } else if (INDEXABLE(Rhs)){
                return Rhs::size;
            } else {
                return 1;
            }
        }

        template<typename Lhs, typename Rhs>
        constexpr bool constexpr_size_eq(){
            if constexpr (INDEXABLE(Lhs) && INDEXABLE(Rhs)){
                return Lhs::size == Rhs::size;
            } else {
                return true;
            }
        }

        template<typename T>
        inline constexpr typename std::enable_if<INDEXABLE(T), typename T::value_type>::type
        expr_at(const T &expr, size_t i) {
            return expr[i];
        }

        template<typename T>
        inline constexpr typename std::enable_if<!INDEXABLE(T), T>::type
        expr_at(const T &expr, size_t i) {
            return expr;
        }

        template<typename Derived>
        struct Expr {
            inline const Derived &cast() const {
                return *static_cast<const Derived *>(this);
            }
        };

        template<typename T>
        std::ostream &operator<<(std::ostream &output, const Expr<T> &var) {
            const T& expr = var.cast();
            const size_t len = expr.size;
            for (size_t i = 0; i < len; ++i) {
                output << expr[i] << " ";
            }
            return output;
        }

        template<typename T>
        std::istream &operator>>(std::istream &input, Expr<T> &var) {
            const T& expr = var.cast();
            const size_t len = expr.size;
            for (size_t i = 0; i < len; ++i) {
                input >> expr[i];
            }
            return input;
        }

        template<typename Operator, typename Unary>
        struct Unary_Expr : public Expr<Unary_Expr<Operator, Unary>> {
        private:
            const Operator &opt_;
            const Unary &unary_;
        public:
            static constexpr size_t size{Unary::size};
            using value_type = typename Unary::value_type;

            Unary_Expr(const Operator &opt, const Unary &unary) : opt_(opt), unary_(unary) {}

            inline constexpr value_type operator[](size_t i) const {
                return opt_(expr_at(unary_, i));
            }
        };

        template<typename Operator, typename Lhs, typename Rhs>
        struct Binary_Expr : public Expr<Binary_Expr<Operator, Lhs, Rhs>> {
            static_assert(constexpr_size_eq<Lhs, Rhs>(), "unequality of array size!");
        private:
            const Operator &opt_;
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            static constexpr size_t size{constexpr_size<Lhs, Rhs>()};
            using value_type = typename Lhs::value_type;

            Binary_Expr(const Operator &opt, const Lhs &lhs, const Rhs &rhs) : opt_(opt), lhs_(lhs), rhs_(rhs) {}

            inline constexpr value_type operator[](size_t i) const {
                return opt_(expr_at(lhs_, i), expr_at(rhs_, i));
            }
        };

#define EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(NAME, FUNC, EXPR)                                                           \
        template<typename UNIQUE_NAME(TMP)>                                                                            \
        const NAME& FUNC(const Expr<UNIQUE_NAME(TMP)> &rhs_expr) {                                                     \
            const UNIQUE_NAME(TMP) &rhs = rhs_expr.cast();                                                             \
            const size_t UNIQUE_NAME(len) = rhs.size;                                                                  \
            for (size_t i = 0; i < UNIQUE_NAME(len); ++i) {                                                            \
                EXPR;                                                                                                  \
             }                                                                                                         \
            return *this;                                                                                              \
        }

#define EXPR_CREATE_UNARY_OPERATION(FUNC, EXPR)                                                                        \
        auto UNIQUE_NAME(OP) = [](const auto& unary ) -> decltype(EXPR) {return (EXPR);};                              \
                                                                                                                       \
        template<typename Unary>                                                                                       \
        inline constexpr Unary_Expr<decltype(UNIQUE_NAME(OP)), Unary>                                                  \
        FUNC(const Expr<Unary> &unary) {                                                                               \
            return Unary_Expr<decltype(UNIQUE_NAME(OP)), Unary>(UNIQUE_NAME(OP), unary.cast());                        \
        }                                                                                                              \

#define EXPR_FILTER(TYPE, ...) typename std::enable_if<!IS_BASE_OF(Expr<TYPE>,TYPE), __VA_ARGS__>::type

#define EXPR_CREATE_BINARY_OPERATION(FUNC, EXPR)                                                                       \
        auto UNIQUE_NAME(OP) =[](const auto &lhs, const auto &rhs)-> decltype(EXPR) {return (EXPR);};                  \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr EXPR_FILTER(Lhs, Binary_Expr<decltype(UNIQUE_NAME(OP)), Lhs, Rhs>)                            \
        FUNC(const Lhs &lhs, const Expr<Rhs> &rhs) {                                                                   \
            return Binary_Expr<decltype(UNIQUE_NAME(OP)), Lhs, Rhs>(UNIQUE_NAME(OP), lhs, rhs.cast());                 \
        }                                                                                                              \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr EXPR_FILTER(Rhs, Binary_Expr<decltype(UNIQUE_NAME(OP)), Lhs, Rhs>)                            \
        FUNC(const Expr<Lhs> &lhs, const Rhs &rhs) {                                                                   \
            return Binary_Expr<decltype(UNIQUE_NAME(OP)), Lhs, Rhs>(UNIQUE_NAME(OP), lhs.cast(), rhs);                 \
        }                                                                                                              \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr Binary_Expr<decltype(UNIQUE_NAME(OP)), Lhs, Rhs>                                              \
        FUNC(const Expr<Lhs> &lhs, const Expr<Rhs> &rhs) {                                                             \
            return Binary_Expr<decltype(UNIQUE_NAME(OP)), Lhs, Rhs>(UNIQUE_NAME(OP), lhs.cast(), rhs.cast());          \
        }                                                                                                              \

        EXPR_CREATE_BINARY_OPERATION(operator+, lhs + rhs);
        EXPR_CREATE_BINARY_OPERATION(operator-, lhs - rhs);
        EXPR_CREATE_BINARY_OPERATION(operator*, lhs * rhs);
        EXPR_CREATE_BINARY_OPERATION(operator/, lhs / rhs);
        EXPR_CREATE_BINARY_OPERATION(operator%, lhs % rhs);
        EXPR_CREATE_BINARY_OPERATION(operator==, lhs == rhs);
        EXPR_CREATE_BINARY_OPERATION(operator!=, lhs != rhs);
        EXPR_CREATE_BINARY_OPERATION(operator<, lhs < rhs);
        EXPR_CREATE_BINARY_OPERATION(operator>, lhs > rhs);
        EXPR_CREATE_BINARY_OPERATION(operator<=, lhs <= rhs);
        EXPR_CREATE_BINARY_OPERATION(operator>=, lhs >= rhs);
        EXPR_CREATE_BINARY_OPERATION(operator&&, lhs && rhs);
        EXPR_CREATE_BINARY_OPERATION(operator||, lhs || rhs);
        EXPR_CREATE_BINARY_OPERATION(operator&, lhs & rhs);
        EXPR_CREATE_BINARY_OPERATION(operator^, lhs ^ rhs);
        EXPR_CREATE_BINARY_OPERATION(operator|, lhs | rhs);
        EXPR_CREATE_BINARY_OPERATION(pow, pow(lhs, rhs));
        EXPR_CREATE_BINARY_OPERATION(exp, exp(lhs, rhs));

        EXPR_CREATE_UNARY_OPERATION(operator+, unary);
        EXPR_CREATE_UNARY_OPERATION(operator-, -unary);
        EXPR_CREATE_UNARY_OPERATION(operator*, *unary);
        EXPR_CREATE_UNARY_OPERATION(operator!, !unary);
        EXPR_CREATE_UNARY_OPERATION(operator~, ~unary);
        EXPR_CREATE_UNARY_OPERATION(abs, unary > 0 ? unary : -unary);
        EXPR_CREATE_UNARY_OPERATION(log, log(unary));
        EXPR_CREATE_UNARY_OPERATION(log10, log10(unary));
        EXPR_CREATE_UNARY_OPERATION(sin, sin(unary));
        EXPR_CREATE_UNARY_OPERATION(cos, cos(unary));
        EXPR_CREATE_UNARY_OPERATION(tan, tan(unary));
        EXPR_CREATE_UNARY_OPERATION(asin, asin(unary));
        EXPR_CREATE_UNARY_OPERATION(acos, acos(unary));
        EXPR_CREATE_UNARY_OPERATION(atan, atan(unary));
        EXPR_CREATE_UNARY_OPERATION(sinh, sinh(unary));
        EXPR_CREATE_UNARY_OPERATION(cosh, cosh(unary));
        EXPR_CREATE_UNARY_OPERATION(tanh, tanh(unary));
    }
}
#endif //SPACEHUB_LAZY_EXPRESSION_H
