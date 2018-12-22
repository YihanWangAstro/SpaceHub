#ifndef SPACEHUB_JOB_SYSTEM_H
#define SPACEHUB_JOB_SYSTEM_H

#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include <deque>

namespace SpaceH {
    namespace MultiThread {



        template<typename T>
        class Opip {
        public:
            explicit Opip(const char *file_name) :
                    file_(file_name, std::fstream::out),
                    thread_(std::thread([&] {
                        while (true) {
                            T data;
                            {
                                std::unique_lock<std::mutex> lock(mutex_);
                                cv_.wait(lock, [&] { return stop_ || !deque_.empty(); });

                                if (stop_)
                                    return;

                                data = std::move(deque_.front());
                                deque_.pop_front();
                            }
                            file_ << data;
                        }
                    })) {}

            ~Opip() {
                {
                    std::lock_guard<std::mutex> lock(mutex_);
                    stop_ = true;
                }
                cv_.notify_all();

                if (thread_.joinable())
                    thread_.join();

                while (!deque_.empty()) {
                    file_ << deque_.front();
                    deque_.pop_front();
                }

            }

            friend void operator<<(Opip &out, T &&tup) {
                {
                    std::lock_guard<std::mutex> lock(out.mutex_);
                    out.deque_.emplace_back(std::forward<T>(tup));
                }
                out.cv_.notify_one();
            }

        private:
            std::deque<T> deque_;
            std::mutex mutex_;
            std::condition_variable cv_;
            std::fstream file_;
            std::thread thread_;
            bool stop_{false};
        };

        template<typename T>
        class Ipip {
        public:
            explicit Ipip(const char *file_name) :
                    file_(file_name, std::fstream::in),
                    thread_(std::thread([&] {
                        while(!file_.eof()){
                            {
                                std::unique_lock<std::mutex> lock(mutex_);
                                cv_load_.wait(lock, [&] { return  stop_ || deque_.empty(); });

                                if(stop_)
                                    break;

                                for(size_t i = 0 ; i < 100 && !file_.eof() ; ++i){
                                    T data;
                                    file_ >> data;
                                    deque_.emplace_back(std::move(data));
                                }
                            }
                            cv_io_.notify_all();
                        }

                        {
                            std::lock_guard<std::mutex> lock(mutex_);
                            eof_ = true;
                        }
                        cv_io_.notify_all();
                    })) {}

            ~Ipip() {
                {
                    std::lock_guard<std::mutex> lock(mutex_);
                    stop_ = true;
                    eof_ = true;
                }
                cv_io_.notify_all();
                cv_load_.notify_one();

                if (thread_.joinable())
                    thread_.join();
            }

            friend bool operator>>(Ipip &in, T &&tup) {
                std::unique_lock<std::mutex> lock(in.mutex_);
                in.cv_io_.wait(lock, [&] { return in.eof_ || !in.deque_.empty(); });

                if (in.eof_ && in.deque_.empty())
                    return false;

                tup = std::move(in.deque_.front());
                in.deque_.pop_front();

                if(in.deque_.empty())
                    in.cv_load_.notify_one();

                return true;
            }

        private:
            std::deque<T> deque_;
            std::mutex mutex_;
            std::condition_variable cv_load_;
            std::condition_variable cv_io_;
            std::fstream file_;
            std::thread thread_;
            bool stop_{false};
            bool eof_{false};
        };

    }
}
#endif //SPACEHUB_JOB_SYSTEM_H
