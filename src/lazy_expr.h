
#ifndef SPACEHUB_LAZY_EXPRESSION_H
#define SPACEHUB_LAZY_EXPRESSION_H

#include "dev_tools.h"
#include "stddef.h"
#include <iostream>

namespace SpaceH {
    namespace Lazy {

        /*
        CREATE_METHOD_CHECK(evaluate);

        template<typename Expr>
        inline constexpr typename std::enable_if<HAS_METHOD(Expr, evaluate, size_t), typename Expr::value_type>::type
        generic_evaluate(const Expr &expr, size_t i) {
            return expr.evaluate(i);
        }

        template<typename Expr>
        inline constexpr typename std::enable_if<!HAS_METHOD(Expr, evaluate, size_t), Expr>::type
        generic_evaluate(const Expr &expr, size_t i) {
            return expr;
        }

        template<typename Uni>
        struct Neg_Expression {
        private:
            const Uni &uni_;
        public:
            using value_type = typename Uni::value_type;

            explicit Neg_Expression(const Uni &uni) : uni_(uni) {}

            inline constexpr value_type evaluate(size_t i) const {
                return -generic_evaluate(uni_, i);
            }
        };

        template<typename Lhs, typename Rhs>
        struct Add_Expression {
        private:
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            using value_type = typename get_value_type<Lhs>::type;

            Add_Expression(const Lhs &lhs, const Rhs &rhs) : lhs_(lhs), rhs_(rhs) {}

            inline constexpr value_type evaluate(size_t i) const {
                return generic_evaluate(lhs_, i) + generic_evaluate(rhs_, i);
            }
        };

        template<typename Lhs, typename Rhs>
        struct Sub_Expression {
        private:
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            using value_type = typename get_value_type<Lhs>::type;

            Sub_Expression(const Lhs &lhs, const Rhs &rhs) : lhs_(lhs), rhs_(rhs) {}

            inline constexpr value_type evaluate(size_t i) const {
                return generic_evaluate(lhs_, i) - generic_evaluate(rhs_, i);
            }
        };

        template<typename Lhs, typename Rhs>
        struct Mul_Expression {
        private:
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            using value_type = typename get_value_type<Lhs>::type;

            Mul_Expression(const Lhs &lhs, const Rhs &rhs) : lhs_(lhs), rhs_(rhs) {}

            inline constexpr value_type evaluate(size_t i) const {
                return generic_evaluate(lhs_, i) * generic_evaluate(rhs_, i);
            }
        };

        template<typename Lhs, typename Rhs>
        struct Div_Expression {
        private:
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            using value_type = typename get_value_type<Lhs>::type;

            Div_Expression(const Lhs &lhs, const Rhs &rhs) : lhs_(lhs), rhs_(rhs) {}

            inline constexpr value_type evaluate(size_t i) const {
                return generic_evaluate(lhs_, i) / generic_evaluate(rhs_, i);
            }
        };

        template<typename Uni>
        inline Neg_Expression<Uni>
        operator-(const Uni &uni) {
            return Neg_Expression<Uni>(uni);
        }

        template<typename Lhs, typename Rhs>
        inline constexpr Add_Expression<Lhs, Rhs>
        operator+(const Lhs &lhs, const Rhs &rhs) {
            return Add_Expression<Lhs, Rhs>(lhs, rhs);
        }

        template<typename Lhs, typename Rhs>
        inline constexpr Sub_Expression<Lhs, Rhs>
        operator-(const Lhs &lhs, const Rhs &rhs) {
            return Sub_Expression<Lhs, Rhs>(lhs, rhs);
        }

        template<typename Lhs, typename Rhs>
        inline constexpr Mul_Expression<Lhs, Rhs>
        operator*(const Lhs &lhs, const Rhs &rhs) {
            return Mul_Expression<Lhs, Rhs>(lhs, rhs);
        }

        template<typename Lhs, typename Rhs>
        inline constexpr Div_Expression<Lhs, Rhs>
        operator/(const Lhs &lhs, const Rhs &rhs) {
            return Div_Expression<Lhs, Rhs>(lhs, rhs);
        }*/

        template<typename Derived>
        struct Expr {
            constexpr Derived &cast() {
                return *static_cast<Derived *>(this);
            }
        };

        template<typename Op, typename Lhs, typename Rhs>
        struct Binary_Expr : public Expr<Binary_Expr<Op, Lhs, Rhs>> {
        private:
            const Lhs &lhs_;
            const Rhs &rhs_;
        public:
            using value_type = typename get_value_type<Lhs>::type;

            Binary_Expr(const Lhs &lhs, const Rhs &rhs) : lhs_(lhs), rhs_(rhs) {}

            inline constexpr value_type evaluate(size_t i) {
                return Op::operation(lhs_.evaluate(i), rhs_.evaluate(i));
            }
        };

#define CREATE_BINARY_OPERATION(NAME, FUN, EXPR)                                                                       \
        template <typename T>                                                                                          \
        struct Op_##NAME {                                                                                             \
            inline constexpr static T operation(T lhs, T rhs){                                                         \
                return EXPR ;                                                                                          \
            }                                                                                                          \
        };                                                                                                             \
                                                                                                                       \
        template<typename Lhs, typename Rhs>                                                                           \
        inline constexpr Binary_Expr<Op_##NAME<typename get_value_type<Lhs>::type>, Lhs, Rhs>                          \
        FUN(const Expr<Lhs> &lhs, const Expr<Rhs> &rhs) {                                                              \
            return Binary_Expr<Op_##NAME<typename get_value_type<Lhs>::type>, Lhs, Rhs>(lhs.cast(), rhs.cast());       \
        }                                                                                                              \


        CREATE_BINARY_OPERATION(Add, operator+, lhs + rhs);

        CREATE_BINARY_OPERATION(Sub, operator-, lhs - rhs);

        CREATE_BINARY_OPERATION(Mul, operator*, lhs * rhs);

        CREATE_BINARY_OPERATION(Div, operator/, lhs / rhs);

        CREATE_BINARY_OPERATION(Pow, pow, pow(lhs, rhs));

        CREATE_BINARY_OPERATION(Exp, exp, exp(lhs, rhs));

        template<typename Op, typename Unary>
        struct Unary_Expr : public Expr<Unary_Expr<Op, Unary>> {
        private:
            const Unary &unary_;
        public:
            using value_type = typename get_value_type<Unary>::type;

            explicit Unary_Expr(const Unary &unary) : unary_(unary) {}

            inline constexpr value_type evaluate(size_t i) {
                return Op::operation(unary_.evaluate(i));
            }
        };

#define CREATE_UNARY_OPERATION(NAME, FUN, EXPR)                                                                        \
        template <typename T>                                                                                          \
        struct Op_##NAME {                                                                                             \
            inline constexpr static T operation(T unary){                                                              \
                return EXPR ;                                                                                          \
            }                                                                                                          \
        };                                                                                                             \
                                                                                                                       \
        template<typename Unary>                                                                                       \
        inline constexpr Unary_Expr<Op_##NAME<typename get_value_type<Unary>::type>, Unary>                            \
        FUN(const Expr<Unary> &unary) {                                                                                \
            return Unary_Expr<Op_##NAME<typename get_value_type<Unary>::type>, Unary>(unary.cast());                   \
        }                                                                                                              \


        CREATE_UNARY_OPERATION(Pos, operator+, unary);

        CREATE_UNARY_OPERATION(Neg, operator-, unary);

        CREATE_UNARY_OPERATION(Abs, abs, unary > 0 ? unary : -unary);

        CREATE_UNARY_OPERATION(Log, log, log(unary));

        CREATE_UNARY_OPERATION(Log10, log10, log10(unary));

        CREATE_UNARY_OPERATION(Sin, sin, sin(unary));

        CREATE_UNARY_OPERATION(Cos, cos, cos(unary));

        CREATE_UNARY_OPERATION(Tan, tan, tan(unary));

        CREATE_UNARY_OPERATION(Asin, asin, asin(unary));

        CREATE_UNARY_OPERATION(Acos, acos, acos(unary));

        CREATE_UNARY_OPERATION(Atan, atan, atan(unary));

        CREATE_UNARY_OPERATION(Sinh, sinh, sinh(unary));

        CREATE_UNARY_OPERATION(Cosh, cosh, cosh(unary));

        CREATE_UNARY_OPERATION(Tanh, tanh, tanh(unary));
    }
}
#endif //SPACEHUB_LAZY_EXPRESSION_H
