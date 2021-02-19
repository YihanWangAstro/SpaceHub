/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file multi-thread.hpp
 *
 * Header file.
 */
#pragma once

#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "../IO.hpp"
#include "../dev-tools.hpp"

/**
 * @namespace space::multi_thread
 * Documentation for Foo here. More docs for Foo here,
 * and down here.
 */
namespace space::multi_thread {

#if defined(_MSC_VER)  // Visual studio
#define thread_local __declspec(thread)
#elif defined(__GCC__)  // GCC
#define thread_local __thread
#endif

    inline const size_t auto_thread =
        (std::thread::hardware_concurrency() > 1) ? std::thread::hardware_concurrency() : 1;

    inline const size_t machine_thread_num =
        (std::thread::hardware_concurrency() > 1) ? std::thread::hardware_concurrency() : 1;

    template <typename Lambda>
    void multi_threads_loop(size_t total_len, size_t thread_num, Lambda &&task) {
        auto len_pth = total_len / thread_num;
        std::vector<std::thread> threads;
        threads.reserve(thread_num);
        for (size_t th_id = 0; th_id < thread_num; ++th_id) {
            size_t begin = th_id * len_pth;
            size_t end = begin + len_pth;
            if (end + len_pth > total_len) end = total_len;
            threads.emplace_back(std::thread(std::forward<Lambda>(task), begin, end));
        }

        for (auto &thread : threads) {
            if (thread.joinable()) thread.join();
        }
    }

    template <typename Callable, typename... Args>
    void multi_thread(size_t thread_num, Callable &&job, Args &&...args) {
        std::vector<std::thread> threads;
        threads.reserve(thread_num);
        for (size_t i = 0; i < thread_num; ++i) {
            threads.emplace_back(std::thread(std::forward<Callable>(job), std::forward<Args>(args)...));
        }

        for (auto &th : threads) {
            if (th.joinable()) th.join();
        }
    }

    template <typename Callable, typename... Args>
    void indexed_multi_thread(size_t thread_num, Callable &&job, Args &&...args) {
        std::vector<std::thread> threads;
        threads.reserve(thread_num);
        for (size_t i = 0; i < thread_num; ++i) {
            threads.emplace_back(std::thread(std::forward<Callable>(job), i, std::forward<Args>(args)...));
        }

        for (auto &th : threads) {
            if (th.joinable()) th.join();
        }
    }

    template <typename... Args>
    void auto_multi_thread(Args &&...args) {
        multi_thread(machine_thread_num, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void auto_indexed_multi_thread(Args &&...args) {
        indexed_multi_thread(machine_thread_num, std::forward<Args>(args)...);
    }

    class ConcurrentFile {
       public:
        inline ConcurrentFile(const char *file_name, std::ios_base::openmode mode);
        inline ConcurrentFile(const std::string &file_name, std::ios_base::openmode mode);

        template <typename Callback, typename... Args>
        auto execute(Callback &&func, Args &&...args);

        inline void flush();

        inline bool eof();

        template <typename U>
        friend ConcurrentFile &operator<<(ConcurrentFile &os, U const &tup);

        template <typename U>
        friend bool operator>>(ConcurrentFile &is, U &tup);

       private:
        std::shared_ptr<std::fstream> file_;
        std::shared_ptr<std::mutex> mutex_;
    };

    ConcurrentFile::ConcurrentFile(const char *file_name, std::ios_base::openmode mode)
        : file_(std::make_shared<std::fstream>(file_name, mode)), mutex_(std::make_shared<std::mutex>()) {
        if (!file_->is_open()) {
            spacehub_abort("Unable to open file: ", file_name);
        }
    }

    ConcurrentFile::ConcurrentFile(const std::string &file_name, std::ios_base::openmode mode)
        : ConcurrentFile(file_name.c_str(), mode) {}

    template <typename Callback, typename... Args>
    auto ConcurrentFile::execute(Callback &&func, Args &&...args) {
        std::lock_guard<std::mutex> lock(*mutex_);
        return func(*file_, std::forward<Args>(args)...);
    }

    void ConcurrentFile::flush() {
        std::lock_guard<std::mutex> lock(*(mutex_));
        file_->flush();
    }

    bool ConcurrentFile::eof() {
        std::lock_guard<std::mutex> lock(*(mutex_));
        return file_->eof();
    }

    // TODO: need to figure out why include IO.hpp doesn't work, thus insert explicitly here.
    template <typename... Args>
    std::ostream &operator<<(std::ostream &out, std::tuple<Args...> const &tup) {
        std::apply(
            [&](auto &&arg, auto &&...args) {
                out << arg;
                (..., (out << ' ' << args));
            },
            tup);
        return out;
    }

    template <typename U>
    ConcurrentFile &operator<<(ConcurrentFile &os, U const &tup) {
        std::lock_guard<std::mutex> lock(*(os.mutex_));
        *(os.file_) << tup;
        return os;
    }

    template <typename U>
    bool operator>>(ConcurrentFile &is, U &tup) {
        std::lock_guard<std::mutex> lock(*(is.mutex_));
        *(is.file_) >> tup;
        bool status = bool(*(is.file_));
        return status;
    }

    inline ConcurrentFile make_thread_safe_fstream(std::string const &name, std::ios_base::openmode mode) {
        return ConcurrentFile(name, mode);
    }
}  // namespace space::multi_thread
