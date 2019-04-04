#ifndef SPACEHUB_SAFE_QUEUE_H
#define SPACEHUB_SAFE_QUEUE_H

#include <atomic>
#include <deque>
#include <mutex>
#include <condition_variable>

namespace SpaceH::MultiThread {

    template<typename T>
    class ConcurrentDeque {
    public:
        ~ConcurrentDeque() {
            stop();
        }

        template<typename ...Args>
        void emplace_back(Args &&... args) {
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

        template<typename ...Args>
        void emplace_front(Args &&... args) {
            {
                std::lock_guard<std::mutex> lock(mutex_);
                deque_.emplace_front(std::forward<Args>(args)...);
            }
            cv_.notify_one();
        }

        T pop_front() noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            while (deque_.empty()) {
                cv_.wait(lock);
            }
            auto data = std::move(deque_.front());
            deque_.pop_front();
            return data;
        }

        T pop_back() noexcept {
            std::lock_guard<std::mutex> lock(mutex_);
            while (deque_.empty()) {
                cv_.wait(lock);
            }
            auto data = std::move(deque_.back());
            deque_.pop_back();
            return data;
        }

    private:
        std::deque<T> deque_;
        std::mutex mutex_;
        std::condition_variable cv_;
        std::atomic_bool stop_{false};
    };
}
#endif //SPACEHUB_SAFE_QUEUE_H
