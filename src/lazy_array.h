
#ifndef SPACEHUB_LAZY_ARRAY_H
#define SPACEHUB_LAZY_ARRAY_H

#include "lazy_expr.h"
namespace SpaceH{
    namespace Lazy{

        constexpr size_t LAZY_ARRAY_DYNAMICAL = 0;

        template<typename Element, size_t Size = LAZY_ARRAY_DYNAMICAL>
        struct Larray {
        public:
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

            inline static constexpr size_t size() {
                return Size;
            }

            inline Element &operator[](size_t i) {
                return data_[i];
            }

            inline const Element &operator[](size_t i) const {
                return data_[i];
            }

            inline Element evaluate(size_t i) const {
                return data_[i];
            }

            template<typename U>
            inline Larray &operator=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] = generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator+=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] += generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator-=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] -= generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator*=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] *= generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator/=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < Size; ++i) {
                    data_[i] /= generic_evaluate(rhs_expr, i);
                }
                return *this;
            }
        private:
            Element data_[Size];
        };


        template<typename Element>
        struct Larray<Element, LAZY_ARRAY_DYNAMICAL> {
        public:
            using value_type = Element;
            Larray() = default;

            template<size_t S>
            Larray(const Element (&src)[S]) : data_(new Element[S],[](Element* p){delete[] p;}), size_(S){
                for (size_t i = 0; i < size_; ++i) {
                    data_.get()[i] = src[i];
                }
            }

            template <typename... T>
            Larray(T ... init_list) : data_{static_cast<Element>(init_list) ...} {
            }

            inline size_t size() {
                return size_;
            }

            inline Element &operator[](size_t i) {
                return data_.get()[i];
            }

            inline const Element &operator[](size_t i) const {
                return data_.get()[i];
            }

            inline Element evaluate(size_t i) const {
                return data_.get()[i];
            }

            template<typename U>
            inline Larray &operator=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < size_; ++i) {
                    data_.get()[i] = generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator+=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < size_; ++i) {
                    data_.get()[i] += generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator-=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < size_; ++i) {
                    data_.get()[i] -= generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator*=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < size_; ++i) {
                    data_.get()[i] *= generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

            template<typename U>
            inline Larray &operator/=(const U &expr) {
                const U &rhs_expr = expr;
                for (size_t i = 0; i < size_; ++i) {
                    data_.get()[i] /= generic_evaluate(rhs_expr, i);
                }
                return *this;
            }

        private:
            std::shared_ptr<Element> data_{nullptr};
            size_t size_{0};
        };

        template<typename Element, size_t Size>
        std::ostream &operator<<(std::ostream &output, const Larray<Element, Size> &larray) {
            size_t size = larray.size();
            for (size_t i = 0; i < size; ++i) {
                output << larray[i] << " ";
            }
            return output;
        }

        template<typename Element, size_t Size>
        std::istream &operator>>(std::istream &input, Larray<Element, Size> &larray) {
            size_t size = larray.size();
            for (size_t i = 0; i < size; ++i) {
                input >> larray[i];
            }
            return input;
        }
    }
}

#endif //SPACEHUB_LAZY_ARRAY_H
