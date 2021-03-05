#pragma once

#include <atomic>
#include <condition_variable>
#include <deque>
#include <mutex>

namespace hub::multi_thread {

    template <typename T>
    class ConcurrentDeque {
       public:
        ~ConcurrentDeque() { stop(); }

        template <typename... Args>
        void emplace_back(Args &&...args) {
            {
                std::lock_guard<std::mutex> lock(mutex_);
                deque_.emplace_back(std::forward<Args>(args)...);
            }
            cv_.notify_one();
        }

        bool empty() {
            std::lock_guard<std::mutex> lock(mutex_);
            return deque_.empty();
        }

        void stop() {
            stop_ = true;
            cv_.notify_all();
        }

        template <typename... Args>
        void emplace_front(Args &&...args) {
            {
                std::lock_guard<std::mutex> lock(mutex_);
                deque_.emplace_front(std::forward<Args>(args)...);
            }
            cv_.notify_one();
        }

        bool try_pop_front(T &data) noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            while (!stop_ && deque_.empty()) {
                cv_.wait(lock);
            }
            if (stop_) return false;

            data = std::move(deque_.front());
            deque_.pop_front();
            return true;
        }

        bool try_pop_back(T &data) noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            while (!stop_ && deque_.empty()) {
                cv_.wait(lock);
            }
            if (stop_) return false;

            data = std::move(deque_.back());
            deque_.pop_back();
            return true;
        }

       private:
        std::deque<T> deque_;
        std::mutex mutex_;
        std::condition_variable cv_;
        std::atomic_bool stop_{false};
    };
}  // namespace hub::multi_thread
