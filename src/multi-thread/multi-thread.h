
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

namespace SpaceH {
    namespace MultiThread {

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

        class LockedFile {
        public:
            LockedFile(const char *file_name, std::ios_base::openmode mode) {
                file_.open(file_name, mode);
                if(!file_.is_open()){
                    std::cout << "Unable to open file: " << std::string(file_name) << "\n";
                    exit(0);
                }
            }

            LockedFile(LockedFile &) = delete;

            template<typename Callback, typename ...Args>
            auto execute(Callback &&func, Args &&...args) {
                std::lock_guard<std::mutex> lock(mutex_);
                return func(file_, std::forward<Args>(args)...);
            }

        private:
            std::fstream file_;
            std::mutex mutex_;
        };

        template<typename T>
        class SharedHolder {
        public:
            SharedHolder(std::shared_ptr<T> shared_obj) : shared_obj_(shared_obj) {}

            template<typename ...Args>
            auto execute(Args &&...args) {
                return shared_obj_->execute(std::forward<Args>(args)...);
            }

        private:
            std::shared_ptr<T> shared_obj_;
        };
    }
}

#endif
