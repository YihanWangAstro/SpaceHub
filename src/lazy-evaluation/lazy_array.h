
#ifndef SPACEHUB_LAZY_ARRAY_H
#define SPACEHUB_LAZY_ARRAY_H

#include "lazy_expr.h"
#include <algorithm>
#include "slice.h"
#include "../multi-thread/multi-thread.hpp"


namespace SpaceH {
    namespace Lazy {

        template<typename T>
        constexpr bool leq_cache_line(size_t len = 1) {
            return sizeof(T) * len <= sizeof(double) * 4;
        }

        template<typename T, size_t Len, bool IsSmall = leq_cache_line<T>(Len)>
        struct Larray : public Expr<Larray<T, Len, IsSmall>> {
        public:
            using value_type = T;

            Larray() : data_(new T[Len]) {}

            Larray(Larray &&src) noexcept : data_(std::move(src.data_)) {
                src.data_ = nullptr;
            }

            template<size_t S>
            Larray(const T (&src)[S]) : data_(new T[Len]) {
                COMPILE_TIME_ASSERT(S <= Len, "Capacity of initializer exceeds array size!");
                std::copy(src, src + S, data_.get());
            }

            template<typename... Args>
            Larray(T arg, Args ... args) : data_(new T[Len]{arg, static_cast<T>(args) ...}) {}

            T const &eval(size_t i) const {
                return data_[i];
            }

            inline static constexpr size_t size() {
                return Len;
            }

            inline const T &operator[](size_t i) const {
                return data_[i];
            }

            inline T &operator[](size_t i) {
                return const_cast<T &>(static_cast<Larray const &>(*this)[i]);
            }

            Larray &operator=(Larray &&src) noexcept {
                data_.reset(std::move(src.data_));
                src.data_ = nullptr;
                return *this;
            }

            Larray &operator=(const T &scalar) {
                std::fill(data_.get(), data_.get() + Len, scalar);
            }

            template<typename U>
            Larray(const Expr<U> &expr) : data_(new T[Len]) {
                const U &src = expr.cast();
                multi_threads_loop(Len, MultiThread::auto_thread, [&](size_t begin, size_t end){
                    for(size_t i = begin; i < end;++i)
                        data_[i] = src.eval(i);
                });
            }

            template<typename U>
            Larray &operator=(const Expr<U> &rhs_expr) {
                const U &rhs = rhs_expr.cast();
                /*auto len_pth = Len / MAX_THREAD_NUM;
                std::vector<std::future<void>> results(MAX_THREAD_NUM);
                for(size_t th_id = 0 ; th_id < MAX_THREAD_NUM; ++th_id){
                    size_t begin = th_id * len_pth;
                    size_t end   = begin + len_pth;
                    if(end + len_pth > Len)
                        end = Len;
                    results[th_id] = POOL.commit([&](size_t begin, size_t end){
                        for(size_t i = begin; i < end; ++i)
                            data_[i] = rhs.eval(i);
                    }, begin, end);
                }

                for(auto & r : results)
                    r.get();*/
                multi_threads_loop(Len, MultiThread::auto_thread, [&](size_t begin, size_t end){
                    for(size_t i = begin; i < end; ++i)
                        data_[i] = rhs.eval(i);
                });
                return *this;
            }

            T max() const {
                T max_value = data_[0];
                for (size_t i = 1; i < Len; ++i) {
                    if (data_[i] > max_value)
                        max_value = data_[i];
                }
                return max_value;
            }

            T min() const {
                T min_value = data_[0];
                for (size_t i = 1; i < Len; ++i) {
                    if (data_[i] < min_value)
                        min_value = data_[i];
                }
                return min_value;
            }

            void sort() {
                std::sort(data_.get(), data_.get() + Len);
            }

            EXPR_CREATE_SLICE_INTERFACE;//create operator[slice(begin,end,stride)]
        private:
            std::unique_ptr<T[]> data_;
        };

        template<typename T, size_t Len>
        struct Larray<T, Len, true> : public Expr<Larray<T, Len, true>> {
        public:
            using value_type = T;

            Larray() = default;

            template<size_t S>
            Larray(const T (&src)[S]) {
                COMPILE_TIME_ASSERT(S <= Len, "Capacity of initializer exceeds array size!");
                std::copy(src, src + S, data_.get());
            }

            template<typename... Args>
            Larray(T arg, Args ... args) : data_{arg, static_cast<T>(args) ...} {}

            T const &eval(size_t i) const {
                return data_[i];
            }

            inline static constexpr size_t size() {
                return Len;
            }

            inline const T &operator[](size_t i) const {
                return data_[i];
            }

            inline T &operator[](size_t i) {
                return const_cast<T &>(static_cast<Larray const &>(*this)[i]);
            }

            Larray &operator=(const T &scalar) {
                std::fill(data_.get(), data_.get() + Len, scalar);
            }

            template<typename U>
            explicit Larray(const Expr<U> &expr) {
                const U &src = expr.cast();
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = src.eval(i);
                }
            }

            template<typename U>
            Larray &operator=(const Expr<U> &rhs_expr) {
                const U &rhs = rhs_expr.cast();
                for (size_t i = 0; i < Len; ++i) {
                    data_[i] = rhs.eval(i);
                }
                return *this;
            }

            T max() const {
                T max_value = data_[0];
                for (size_t i = 1; i < Len; ++i) {
                    if (data_[i] > max_value)
                        max_value = data_[i];
                }
                return max_value;
            }

            T min() const {
                T min_value = data_[0];
                for (size_t i = 1; i < Len; ++i) {
                    if (data_[i] < min_value)
                        min_value = data_[i];
                }
                return min_value;
            }

            void sort() {
                std::sort(data_.get(), data_.get() + Len);
            }

            EXPR_CREATE_SLICE_INTERFACE;//create operator[slice(begin,end,stride)]
        private:
            T data_[Len];
        };


    }
}

#endif //SPACEHUB_LAZY_ARRAY_H
