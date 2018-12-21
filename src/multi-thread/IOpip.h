#ifndef SPACEHUB_JOB_SYSTEM_H
#define SPACEHUB_JOB_SYSTEM_H

#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include "concurrent-queue.h"

namespace SpaceH {
    namespace MultiThread {

        template<typename T>
        class Opip {
        public:
            explicit Opip(const char *file_name) :
                file_(file_name, std::fstream::out),
                thread_(std::thread([&]{
                    T data;
                    while(pip_.pop_front(data)){
                        file_ << data;
                    }
                })){}

            ~Opip(){
                stop_ = true;
                if (thread_.joinable())
                    thread_.join();

                while(!pip_.empty()){
                    file_ << pip_.pop_front();
                }
            }
            friend Opip& operator<<(Opip& out, T&& tup){
                out.pip_.emplace_back(std::forward<T>(tup));
            }
        private:
            ConcurrentDeque<T> pip_;
            std::fstream file_;
            std::atomic_bool stop_{false};
            std::thread thread_;
        };


        /*template<typename T>
        class IOpip {
        public:
            IOpip(const char *file_name, std::ios_base::openmode mode) :
                file_(std::make_shared<std::fstream>(file_name, mode)),
                pip_(std::make_shared<ConcurrentDeque<T>>()),
                stop_(std::make_shared<std::atomic<bool>>(false)){}

            void stop(){
                *stop_ = true;
            }
        private:
            std::shared_ptr<ConcurrentDeque<T>> pip_;
            std::shared_ptr<std::fstream> file_;
            std::shared_ptr<std::atomic<bool>> stop_;
        };

        template <typename T>
        class Ipip : public IOpip<T>{
        public:
            explicit Ipip(const char *file_name) : IOpip<T>(file_name, std::fstream::in){}

            explicit Ipip(const std::string& file_name) : Ipip(file_name.c_str()){}

            void read_one_from_file(){
                T data;
                *(this->file_) >> data;
                this->pip_->emplace_bak(std::move(data));
            }

            friend Ipip& operator>>(Ipip& in, T&& tup){
                tup = std::move(in.pip_->pop_front());
            }
        };

        template <typename T>
        class Opip : public IOpip<T>{
        public:
            explicit Opip(const char *file_name) : IOpip<T>(file_name, std::fstream::in){}

            explicit Opip(const std::string& file_name) : Opip(file_name.c_str()){}

            void start(){
                while(!(this->stop_)){
                    *(this->file_) << this->pip_->pop_front();
                }
            }

            void write_one_to_file(){
                *(this->file_) << this->pip_->pop_front();
            }

            friend Opip& operator<<(Opip& out, T&& tup){
                out.pip_->emplace_back(std::forward<T>(tup));
            }
        };*/
    }
}
#endif //SPACEHUB_JOB_SYSTEM_H
