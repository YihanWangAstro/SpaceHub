#ifndef SPACEHUB_JOB_SYSTEM_H
#define SPACEHUB_JOB_SYSTEM_H

#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include <queue>
#include "dev-tools.h"

namespace SpaceH::MultiThread {
    template<typename T>
    class Opip {
    public:
        explicit Opip(const char *file_name) :
                file_(file_name, std::fstream::out),
                thread_(std::thread([&] {
                    while (true) {
                        std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
                        lock.lock();
                        cv_.wait(lock, [&] { return stop_ || !queue_.empty(); });

                        if (stop_)
                            return;

                        T data = std::move(queue_.front());
                        queue_.pop_front();
                        lock.unlock();

                        file_ << data;
                    }

                    while (!queue_.empty()) {
                        file_ << queue_.front();
                        queue_.pop_front();
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
        }

        friend void operator<<(Opip &out, T &&tup) {
            {
                std::lock_guard<std::mutex> lock(out.mutex_);
                out.queue_.emplace_back(std::forward<T>(tup));
            }
            out.cv_.notify_one();
        }

    private:
        std::queue<T> queue_;
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
                    while (!file_.eof()) {
                        {
                            std::unique_lock<std::mutex> lock(mutex_);
                            cv_load_.wait(lock, [&] { return stop_ || queue_.empty(); });

                            if (stop_)
                                break;

                            for (size_t i = 0; i < 100 && !file_.eof(); ++i) {
                                T data;
                                file_ >> data;
                                queue_.emplace_back(std::move(data));
                            }
                        }
                        cv_io_.notify_all();
                    }

                    {
                        std::lock_guard<std::mutex> lock(mutex_);
                        eoq_ = true;
                    }
                    cv_io_.notify_all();
                })) {}

        ~Ipip() {
            {
                std::lock_guard<std::mutex> lock(mutex_);
                stop_ = true;
                eoq_ = true;
            }
            cv_io_.notify_all();
            cv_load_.notify_one();

            if (thread_.joinable())
                thread_.join();
        }

        friend bool operator>>(Ipip &in, T &&tup) {
            std::unique_lock<std::mutex> lock(in.mutex_);
            in.cv_io_.wait(lock, [&] { return in.eoq_ || !in.queue_.empty(); });

            if (in.eoq_ && in.queue_.empty())
                return false;

            tup = std::move(in.queue_.front());
            in.queue_.pop_front();

            if (in.queue_.empty())
                in.cv_load_.notify_one();

            return true;
        }

    private:
        std::queue<T> queue_;
        std::mutex mutex_;
        std::condition_variable cv_load_;
        std::condition_variable cv_io_;
        std::fstream file_;
        std::thread thread_;
        bool stop_{false};
        bool eoq_{false};
    };

    template<typename T>
    class COstream {
    public:
        explicit COstream(const char *file_name) : pip_ptr_(std::make_unique<Opip<T>>(file_name)) {}

        explicit COstream(std::string &file_name) : COstream(file_name.c_str()) {}

        COstream(COstream const &) = delete;

        COstream &operator=(COstream const &) = delete;

        COstream(COstream &&other) noexcept : pip_ptr_(std::move(other.pip_ptr_)) {}

        COstream &operator=(COstream &&other) noexcept {
            pip_ptr_ = std::move(other.pip_ptr_);
            return *this;
        }

        friend void operator<<(COstream &out, T const &tup) {
            *(out.pip_ptr_) << tup;
        }

    private:
        std::unique_ptr<Opip<T>> pip_ptr_;
    };

    template<typename T>
    class CIstream {
    public:
        explicit CIstream(const char *file_name) : pip_ptr_(std::make_shared<Ipip<T>>(file_name)) {}

        explicit CIstream(std::string &file_name) : CIstream(file_name.c_str()) {}

        friend void operator>>(CIstream &in, T &&tup) {
            *(in.pip_ptr_) >> std::forward<T>(tup);
        }

    private:
        std::shared_ptr<Ipip<T>> pip_ptr_;
    };

}
#endif //SPACEHUB_JOB_SYSTEM_H
