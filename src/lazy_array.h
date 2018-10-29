
#ifndef SPACEHUB_LAZY_ARRAY_H
#define SPACEHUB_LAZY_ARRAY_H

#include "lazy_expr.h"
namespace SpaceH{
    namespace Lazy{
/*
        template<typename Element, size_t Size>
        struct Larray : public Expr<Larray<Element, Size>>{
        public:
            static constexpr size_t size{Size};
            using value_type = Element;

            Larray() = default;

            template<size_t S>
            Larray(const Element (&src)[S]) {
                COMPILE_TIME_ASSERT(S == Size, "Size inconsistency of initializer and array!");
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] = src[i];
                }
            }
            template <typename... T>
            Larray(T ... init_list) : data_{static_cast<Element>(init_list) ...} {
            }

            inline static constexpr size_t len() {
                return Size;
            }

            inline Element &operator[](size_t i) {
                return data_[i];
            }

            inline const Element &operator[](size_t i) const {
                return data_[i];
            }

            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator=,  data_[i]  = rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator+=, data_[i] += rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator-=, data_[i] -= rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator*=, data_[i] *= rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator/=, data_[i] /= rhs[i]);

        private:
            Element data_[Size];
        };
        */

        template<typename Element, size_t Size>
        struct Larray : public Expr<Larray<Element, Size>>{
        public:
            static constexpr size_t size{Size};
            using value_type = Element;

            Larray() : data_(new Element[Size]){};

            template<size_t S>
            Larray(const Element (&src)[S]) : data_(new Element[Size]) {
                COMPILE_TIME_ASSERT(S == Size, "Size inconsistency of initializer and array!");
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] = src[i];
                }
            }
            template <typename... T>
            Larray(T ... init_list) : data_(new Element[Size]{static_cast<Element>(init_list) ...}) {}

            inline static constexpr size_t len() {
                return Size;
            }

            inline Element &operator[](size_t i) {
                return data_[i];
            }

            inline const Element &operator[](size_t i) const {
                return data_[i];
            }

            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator=,  data_[i]  = rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator+=, data_[i] += rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator-=, data_[i] -= rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator*=, data_[i] *= rhs[i]);
            EXPR_CREATE_LINEAR_ASSIGN_OPERATOR(Larray, operator/=, data_[i] /= rhs[i]);

        private:
            std::unique_ptr<Element[]> data_;
        };
    }
}

#endif //SPACEHUB_LAZY_ARRAY_H
