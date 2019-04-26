#ifndef SPACEHUB_MKL_LINKED_ARRAY_H
#define SPACEHUB_MKL_LINKED_ARRAY_H

#include "lazy_expr.h"
#include <string.h>
#include "/opt/intel/mkl/include/mkl.h"
#include "../tools/timer.hpp"

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
            return *this;
        }

        MKLArray(MKLArray &&src) noexcept {
            data_ = src.data_;
            src.data_ = nullptr;
        }

        MKLArray &operator=(MKLArray &&src) noexcept {
            delete[] data_;
            data_ = src.data_;
            src.data_ = nullptr;
            return *this;
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

}

#endif //SPACEHUB_MKL_LINKED_ARRAY_H
