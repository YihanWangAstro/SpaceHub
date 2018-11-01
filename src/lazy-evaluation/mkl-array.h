#ifndef SPACEHUB_MKL_LINKED_ARRAY_H
#define SPACEHUB_MKL_LINKED_ARRAY_H

#include "lazy_expr.h"
#include <string.h>
#include "/opt/intel/mkl/include/mkl.h"
#include "../tools/timer.h"
namespace SpaceH {

    template<typename T, size_t Size>
    struct MKLArray {
        static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                      "Only float and double are allowed.");

        MKLArray() : data_(new T[Size]) {}

        ~MKLArray() { delete[] data_; }

        MKLArray(const MKLArray &src) : data_(new T[Size]) {
            memcpy(data_, src.data_, sizeof(T) * Size);
        }

        MKLArray &operator=(MKLArray const &src) {
            memcpy(data_, src.data_, sizeof(T) * Size);
            std::cout << "cp\n";
            return *this;
        }

        MKLArray(MKLArray &&src) {
            data_ = src.data_;
            src.data_ = nullptr;
            std::cout << "mv ctor\n";
        }

        MKLArray &operator=(MKLArray &&src) {
            delete[] data_;
            data_ = src.data_;
            src.data_ = nullptr;
            //std::cout << "mv\n";
        }

        T const &operator[](size_t i) const {
            return data_[i];
        }

        T &operator[](size_t i) {
            return const_cast<T &>(static_cast<MKLArray const &>(*this)[i]);
        }

        T *get() const {
            return data_;
        }

    private:
        T *data_{nullptr};
    };

#define MKL_CREATE_BINARY_OPERATION(FUNC, INTERFACE, TYPE)                                                             \
    template <size_t Size>                                                                                             \
    inline MKLArray<TYPE,Size> FUNC(MKLArray<TYPE,Size> const &lhs, MKLArray<TYPE,Size> const &rhs){                   \
        MKLArray<TYPE,Size> dst;                                                                                       \
        INTERFACE(static_cast<MKL_INT>(Size), lhs.get(), rhs.get(), dst.get());                                        \
        return dst;                                                                                                    \
    }

    MKL_CREATE_BINARY_OPERATION(operator+, vdAdd, double);

    MKL_CREATE_BINARY_OPERATION(operator-, vdSub, double);

    MKL_CREATE_BINARY_OPERATION(operator*, vdMul, double);

    MKL_CREATE_BINARY_OPERATION(operator/, vdDiv, double);

    MKL_CREATE_BINARY_OPERATION(operator+, vsAdd, float);

    MKL_CREATE_BINARY_OPERATION(operator-, vsSub, float);

    MKL_CREATE_BINARY_OPERATION(operator*, vsMul, float);

    MKL_CREATE_BINARY_OPERATION(operator/, vsDiv, float);

    namespace Lazy {

        template<typename T, size_t Size>
        struct LazyMKL : public Expr<LazyMKL< T, Size>> {
        public:
            static constexpr size_t len{1};
            using value_type = MKLArray<T, Size>;

            LazyMKL() : mkl_array_ptr_(new MKLArray<T, Size>){}
            ~LazyMKL(){delete mkl_array_ptr_;}
            LazyMKL(const LazyMKL& src) : mkl_array_ptr_(new MKLArray<T, Size>){
                *mkl_array_ptr_ = *src.mkl_array_ptr_;
            }
            LazyMKL& operator=(const LazyMKL& src){
                *mkl_array_ptr_ = *src.mkl_array_ptr_;
                std::cout << "wr =\n";
                return *this;
            }
            LazyMKL(LazyMKL&& src) {
                mkl_array_ptr_ = src.mkl_array_ptr_;
                src.mkl_array_ptr_ = nullptr;
            }
            LazyMKL& operator=(LazyMKL&& src){
                delete mkl_array_ptr_;
                mkl_array_ptr_ = src.mkl_array_ptr_;
                src.mkl_array_ptr_ = nullptr;
                std::cout << "wr mv=\n";
                return *this;
            }

            value_type const &eval(size_t i) const {
                return *mkl_array_ptr_;
            }

            value_type &dst(size_t i) {
                return *mkl_array_ptr_;
            }

            T const &operator[](size_t i) const {
                return (*mkl_array_ptr_)[i];
            }

            T &operator[](size_t i) {
                return const_cast<T &>(static_cast<LazyMKL const &>(*this)[i]);
            }

            EXPR_CREATE_LINEAR_ASSIGN_OPERATORS;//create operator=,+=,-=,*=,/= for Expr<T> as rhs.
        private:
            MKLArray<T, Size>* mkl_array_ptr_;
        };
    }
}

#endif //SPACEHUB_MKL_LINKED_ARRAY_H
