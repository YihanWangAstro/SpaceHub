#pragma once

#include <cmath>

#include "../dev-tools.hpp"

namespace space::lazy {

#define IS_EXPR(TYPE) IS_BASE_OF(Expr<TYPE>, TYPE)

    template <typename Derived>
    struct Expr {
        inline const Derived &cast() const { return *static_cast<const Derived *>(this); }
    };

    template <bool IsExpr, typename T, typename... Args>
    inline const auto eval_at(std::enable_if_t<IsExpr, T> const &expr, Args &&...args) {
        return expr.eval(std::forward<Args>(args)...);
    }

    template <bool IsExpr, typename T, typename... Args>
    inline const auto eval_at(std::enable_if_t<!IsExpr, T> const &expr, Args &&...args) {
        return expr;
    }

    template <typename Operator, typename Unary>
    struct Unary_Expr : public Expr<Unary_Expr<Operator, Unary>> {
       private:
        const Operator &opt_;
        const Unary &unary_;

       public:
        Unary_Expr(const Operator &opt, const Unary &unary) : opt_(opt), unary_(unary) {}

        template <typename... Args>
        inline auto eval(Args &&...args) const {
            return opt_(eval_at<IS_EXPR(Unary), Unary, Args...>(unary_, std::forward<Args>(args)...));
        }
    };

    template <typename Operator, typename Lhs, typename Rhs>
    struct Binary_Expr : public Expr<Binary_Expr<Operator, Lhs, Rhs>> {
       private:
        const Operator &opt_;
        const Lhs &lhs_;
        const Rhs &rhs_;

       public:
        Binary_Expr(const Operator &opt, const Lhs &lhs, const Rhs &rhs) : opt_(opt), lhs_(lhs), rhs_(rhs) {}

        template <typename... Args>
        inline auto eval(Args &&...args) const {
            return opt_(eval_at<IS_EXPR(Lhs), Lhs, Args...>(lhs_, std::forward<Args>(args)...),
                        eval_at<IS_EXPR(Rhs), Rhs, Args...>(rhs_, std::forward<Args>(args)...));
        }
    };

#define EXPR_CREATE_UNARY_OPERATION(FUNC, EXPR)                                             \
    inline auto UNIQ(OP) = [](const auto &unary) -> decltype(EXPR) { return (EXPR); };      \
                                                                                            \
    template <typename Unary>                                                               \
    inline constexpr Unary_Expr<decltype(UNIQ(OP)), Unary> FUNC(const Expr<Unary> &unary) { \
        return Unary_Expr<decltype(UNIQ(OP)), Unary>(UNIQ(OP), unary.cast());               \
    }

#define EXPR_FILTER(TYPE, ...) typename std::enable_if_t<!IS_EXPR(TYPE), ##__VA_ARGS__>

#define EXPR_CREATE_BINARY_OPERATION(FUNC, EXPR)                                                                  \
    inline auto UNIQ(OP) = [](const auto &lhs, const auto &rhs) -> decltype(EXPR) { return (EXPR); };             \
                                                                                                                  \
    template <typename Lhs, typename Rhs>                                                                         \
    inline constexpr EXPR_FILTER(Lhs, Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>)                                  \
        FUNC(const Lhs &lhs, const Expr<Rhs> &rhs) {                                                              \
        return Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>(UNIQ(OP), lhs, rhs.cast());                              \
    }                                                                                                             \
                                                                                                                  \
    template <typename Lhs, typename Rhs>                                                                         \
    inline constexpr EXPR_FILTER(Rhs, Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>)                                  \
        FUNC(const Expr<Lhs> &lhs, const Rhs &rhs) {                                                              \
        return Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>(UNIQ(OP), lhs.cast(), rhs);                              \
    }                                                                                                             \
                                                                                                                  \
    template <typename Lhs, typename Rhs>                                                                         \
    inline constexpr Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs> FUNC(const Expr<Lhs> &lhs, const Expr<Rhs> &rhs) { \
        return Binary_Expr<decltype(UNIQ(OP)), Lhs, Rhs>(UNIQ(OP), lhs.cast(), rhs.cast());                       \
    }

    EXPR_CREATE_BINARY_OPERATION(operator+, lhs + rhs);
    EXPR_CREATE_BINARY_OPERATION(operator-, lhs - rhs);
    EXPR_CREATE_BINARY_OPERATION(operator*, lhs *rhs);
    EXPR_CREATE_BINARY_OPERATION(operator/, lhs / rhs);
    EXPR_CREATE_BINARY_OPERATION(operator%, lhs % rhs);
    EXPR_CREATE_BINARY_OPERATION(operator==, lhs == rhs);
    EXPR_CREATE_BINARY_OPERATION(operator!=, lhs != rhs);
    EXPR_CREATE_BINARY_OPERATION(operator<, lhs < rhs);
    EXPR_CREATE_BINARY_OPERATION(operator>, lhs > rhs);
    EXPR_CREATE_BINARY_OPERATION(operator<=, lhs <= rhs);
    EXPR_CREATE_BINARY_OPERATION(operator>=, lhs >= rhs);
    EXPR_CREATE_BINARY_OPERATION(operator&&, lhs &&rhs);
    EXPR_CREATE_BINARY_OPERATION(operator||, lhs || rhs);
    EXPR_CREATE_BINARY_OPERATION(operator&, lhs &rhs);
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
}  // namespace space::lazy
