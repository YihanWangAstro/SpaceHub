
#ifndef SPACEHUB_LAZY_EXPRESSION_H
#define SPACEHUB_LAZY_EXPRESSION_H

#include "../dev_tools.h"
#include "stddef.h"
#include <math.h>

namespace SpaceH {
    namespace Lazy {

#define IS_EXPR(TYPE) IS_BASE_OF(Expr<TYPE>,TYPE)

        template<typename Derived>
        struct Expr {
            inline const Derived &cast() const {
                return *static_cast<const Derived *>(this);
            }
        };

        template<typename Lhs, typename Rhs>
        constexpr size_t constexpr_len(){
            if (IS_EXPR(Lhs)){
                return Lhs::len;
            } else if (IS_EXPR(Rhs)){
                return Rhs::len;
            } else {
                return 1;
            }
        }

        template<typename Lhs, typename Rhs>
        constexpr bool constexpr_len_eq(){
            if (IS_EXPR(Lhs) && IS_EXPR(Rhs)){
                return Lhs::len == Rhs::len;
            } else {
                return true;
            }
        }

        template<typename T>
        inline const typename std::enable_if<IS_EXPR(T), typename T::value_type>::type
        eval_at(const T &expr, size_t i) {
            return expr.eval(i);
        }

        template<typename T>
        inline const typename std::enable_if<!IS_EXPR(T), T>::type
        eval_at(const T &expr, size_t i) {
            return expr;
        }



        template<typename T>
        std::ostream &operator<<(std::ostream &output, const Expr<T> &var) {
            const T& expr = var.cast();
            const size_t len = expr.len;
            for (size_t i = 0; i < len; ++i) {
                output << expr.eval(i) << " ";
            }
            return output;
        }

        template<typename T>
        std::istream &operator>>(std::istream &input, Expr<T> &var) {
            const T& expr = var.cast();
            const size_t len = expr.len;
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
            static constexpr size_t len{Unary::len};
            using value_type = typename Unary::value_type;

            Unary_Expr(const Operator &opt, const Unary &unary) : opt_(opt), unary_(unary) {}

            inline const value_type eval(size_t i) const {
                return opt_(eval_at(unary_, i));
            }
        };

        template<typename Operator, typename Lhs, typename Rhs>
        struct Binary_Expr : public Expr<Binary_Expr<Operator, Lhs, Rhs>> {
            static_assert(constexpr_len_eq<Lhs, Rhs>(), "unequality of array size!");
        private:
            const Operator &opt_;
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            static constexpr size_t len{constexpr_len<Lhs, Rhs>()};
            using value_type = typename Lhs::value_type;

            Binary_Expr(const Operator &opt, const Lhs &lhs, const Rhs &rhs) : opt_(opt), lhs_(lhs), rhs_(rhs) {}

            inline const value_type eval(size_t i) const {
                return opt_(eval_at(lhs_, i), eval_at(rhs_, i));
            }
        };

#define EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(FUNC, EXPR)                                                                 \
        template<typename UNIQ(TMP)>                                                                                   \
        auto FUNC(const Expr<UNIQ(TMP)> &rhs_expr)-> REF_TYPE_OF_SELF {                                                \
            const UNIQ(TMP) &rhs = rhs_expr.cast();                                                                    \
            const size_t UNIQ(len) = rhs.len;                                                                          \
            for (size_t i = 0; i < UNIQ(len); ++i) {                                                                   \
                EXPR;                                                                                                  \
             }                                                                                                         \
            return *this;                                                                                              \
        }

#define EXPR_CREATE_UNARY_OPERATION(FUNC, EXPR)                                                                        \
        auto UNIQ(OP) = [](const auto& unary ) -> decltype(EXPR) {return (EXPR);};                                     \
                                                                                                                       \
        template<typename Unary>                                                                                       \
        inline constexpr Unary_Expr<decltype(UNIQ(OP)), Unary>                                                         \
        FUNC(const Expr<Unary> &unary) {                                                                               \
            return Unary_Expr<decltype(UNIQ(OP)), Unary>(UNIQ(OP), unary.cast());                                      \
        }                                                                                                              \

#define EXPR_FILTER(TYPE, ...) typename std::enable_if<!IS_EXPR(TYPE), __VA_ARGS__>::type

#define EXPR_CREATE_BINARY_OPERATION(FUNC, EXPR)                                                                       \
        auto UNIQ(OP) =[](const auto &lhs, const auto &rhs)-> decltype(EXPR) {return (EXPR);};                         \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr EXPR_FILTER(Lhs, Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>)                                   \
        FUNC(const Lhs &lhs, const Expr<Rhs> &rhs) {                                                                   \
            return Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>(UNIQ(OP), lhs, rhs.cast());                               \
        }                                                                                                              \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr EXPR_FILTER(Rhs, Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>)                                   \
        FUNC(const Expr<Lhs> &lhs, const Rhs &rhs) {                                                                   \
            return Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>(UNIQ(OP), lhs.cast(), rhs);                               \
        }                                                                                                              \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>                                                     \
        FUNC(const Expr<Lhs> &lhs, const Expr<Rhs> &rhs) {                                                             \
            return Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>(UNIQ(OP), lhs.cast(), rhs.cast());                        \
        }                                                                                                              \

#define EXPR_CREATE_LINEAR_ASSIGN_OPERATORS                                                                            \
        EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(operator=,  (*this).dst(i)  = rhs.eval(i));                                 \
        EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(operator+=, (*this).dst(i) += rhs.eval(i));                                 \
        EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(operator-=, (*this).dst(i) -= rhs.eval(i));                                 \
        EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(operator*=, (*this).dst(i) *= rhs.eval(i));                                 \
        EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(operator/=, (*this).dst(i) /= rhs.eval(i));                                 \

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
