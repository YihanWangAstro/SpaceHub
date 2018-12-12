
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
            
            ~LockedFile(){
                file_.close();
            }

            LockedFile(LockedFile &) = delete;

            template<typename Callback, typename ...Args>
            auto execute(Callback &&func, Args &&...args) {
                std::lock_guard<std::mutex> lock(mutex_);
                return func(file_, std::forward<Args>(args)...);
            }

            template<typename ...Args>
            friend LockedFile& operator<<(LockedFile& os, std::tuple<Args...>&& tup){
                std::lock_guard<std::mutex> lock(os.mutex_);
                os.print(std::forward<decltype(tup)>(tup), std::make_index_sequence<sizeof...(Args)>());
                return os;
            }
            
            template<typename ...Args>
            friend LockedFile& operator>>(LockedFile& is, std::tuple<Args...>&& tup){
                std::lock_guard<std::mutex> lock(is.mutex_);
                is.input(std::forward<decltype(tup)>(tup), std::make_index_sequence<sizeof...(Args)>());
                return is;
            }

        private:
            template<class Tup, size_t... I>
            void print(const Tup& tup, std::index_sequence<I...>)
            {
                (..., (file_ << std::get<I>(tup)));
            }

            template<class Tup, size_t... I>
            void input(const Tup& tup, std::index_sequence<I...>)
            {
                (..., (file_ >> std::get<I>(tup)));
            }
        private:
            std::fstream file_;
            std::mutex mutex_;
        };

#define PACK(...) std::forward_as_tuple(__VA_ARGS__)
            
        template<typename T>
        class SharedHolder {
        public:
            SharedHolder(std::shared_ptr<T> shared_obj) : shared_obj_(shared_obj) {}

            template<typename ...Args>
            auto execute(Args &&...args) {
                return shared_obj_->execute(std::forward<Args>(args)...);
            }
            
            template<typename ...Args>
            friend SharedHolder& operator<<(SharedHolder& os , std::tuple<Args...>&& tup){
                *(os.shared_obj_) << std::forward<decltype(tup)>(tup);
                return os;
            }

            template<typename ...Args>
            friend SharedHolder& operator>>(SharedHolder& is , std::tuple<Args...>&& tup){
                *(is.shared_obj_) >> std::forward<decltype(tup)>(tup);
                return is;
            }
        private:
            
        private:
            std::shared_ptr<T> shared_obj_;
        };

        using ThreadSharedFile = SharedHolder<LockedFile>;

        ThreadSharedFile make_thread_safe_fstream(const char *file_name, std::ios_base::openmode mode){
            return ThreadSharedFile(std::make_shared<LockedFile>(file_name,mode));
        }
    }
}

#endif
