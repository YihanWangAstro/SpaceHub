
#ifndef SPACEHUB_SYN_FILE_H
#define SPACEHUB_SYN_FILE_H

#include <fstream>
#include <memory>
#include <string>
#include <stdio.h>
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>
#include <tuple>
#include "dev-tools.h"

namespace space::multiThread {

    const size_t auto_thread = (std::thread::hardware_concurrency() > 1) ? std::thread::hardware_concurrency() : 1;

    template<typename Lambda>
    void multi_threads_loop(size_t total_len, size_t thread_num, Lambda &&task) {
        auto len_pth = total_len / thread_num;
        std::vector<std::thread> threads;
        threads.reserve(thread_num);
        for (size_t th_id = 0; th_id < thread_num; ++th_id) {
            size_t begin = th_id * len_pth;
            size_t end = begin + len_pth;
            if (end + len_pth > total_len)
                end = total_len;
            threads.emplace_back(std::thread(std::forward<Lambda>(task), begin, end));
        }

        for (auto &thread : threads) {
            if (thread.joinable())
                thread.join();
        }
    }

    template<typename Callable, typename ...Args>
    void multi_thread_run(size_t thread_num, Callable &&job, Args &&...args) {
        std::vector<std::thread> threads;
        threads.reserve(thread_num);
        for (size_t i = 0; i < thread_num; ++i) {
            threads.emplace_back(std::thread(std::forward<Callable>(job), std::forward<Args>(args)...));
        }

        for (auto &th : threads) {
            if (th.joinable())
                th.join();
        }
    }

    template<typename ...Args>
    void auto_multi_thread(Args &&... args) {
        multi_thread_run(auto_thread, std::forward<Args>(args)...);
    }

    class ConcurrentFile {
    public:
        ConcurrentFile(const char *file_name, std::ios_base::openmode mode) :
                file_(std::make_shared<std::fstream>(file_name, mode)), mutex_(std::make_shared<std::mutex>()) {
            if (!file_->is_open()) {
                SPACEHUB_ABORT("Unable to open file: ", file_name);
            }
        }

        ConcurrentFile(const std::string &file_name, std::ios_base::openmode mode) :
                ConcurrentFile(file_name.c_str(), mode) {}

        template<typename Callback, typename ...Args>
        auto execute(Callback &&func, Args &&...args) {
            std::lock_guard<std::mutex> lock(*mutex_);
            return func(*file_, std::forward<Args>(args)...);
        }

        template<typename U>
        friend ConcurrentFile &operator<<(ConcurrentFile &os, U &&tup) {
            std::lock_guard<std::mutex> lock(*(os.mutex_));
            *(os.file_) << std::forward<U>(tup);
            return os;
        }

        template<typename U>
        friend ConcurrentFile &operator>>(ConcurrentFile &is, U &&tup) {
            std::lock_guard<std::mutex> lock(*(is.mutex_));
            *(is.file_) >> std::forward<U>(tup);
            return is;
        }

    private:
        std::shared_ptr<std::fstream> file_;
        std::shared_ptr<std::mutex> mutex_;
    };

    ConcurrentFile make_thread_safe_fstream(std::string const& name, std::ios_base::openmode mode){
        return ConcurrentFile(name, mode);
    }
}

#endif
