#pragma once

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

namespace space::multi_thread {
    /*---------------------------------------------------------------------------*\
         Class ThreadPool Declaration
    \*---------------------------------------------------------------------------*/
    class ThreadPool {
       public:
        using Task = std::function<void()>;

        template <typename Func, typename... Args>
        auto commit(Func &&func, Args &&...args);

        static ThreadPool &get_instance() {
            static ThreadPool instance{
                (std::thread::hardware_concurrency() > 1) ? std::thread::hardware_concurrency() : 1, 1000};
            return instance;
        }

       private:
        ThreadPool() = delete;

        ThreadPool(ThreadPool const &) = delete;

        ThreadPool &operator=(ThreadPool const &) = delete;

        ThreadPool(size_t max_thread, size_t max_tasks);

        ~ThreadPool();

        std::vector<std::thread> workers_;
        std::queue<Task> tasks_;
        std::mutex lock_;
        std::condition_variable cv_;
        std::atomic_bool stop_;
        const size_t max_tasks_;
    };

    /*---------------------------------------------------------------------------*\
         Class ThreadPool Implementation
    \*---------------------------------------------------------------------------*/
    ThreadPool::ThreadPool(size_t max_thread, size_t max_tasks) : stop_(false), max_tasks_(max_tasks) {
        for (size_t i = 0; i < max_thread; ++i) {
            workers_.emplace_back([&]() {
                Task task;
                while (true) {
                    {
                        std::unique_lock<std::mutex> lock(lock_);
                        while (!stop_ && tasks_.empty()) {
                            cv_.wait(lock);
                        }
                        if (stop_) return;
                        task = tasks_.front();
                        tasks_.pop();
                    }
                    task();
                }
            });
        }
    }

    ThreadPool::~ThreadPool() {
        stop_ = true;
        cv_.notify_all();
        for (auto &worker : workers_) {
            if (worker.joinable()) worker.join();
        }
    }

    template <typename Func, typename... Args>
    auto ThreadPool::commit(Func &&func, Args &&...args) {
        using ReturnType = decltype(func(args...));
        using Package = std::packaged_task<ReturnType()>;

        auto job = std::make_shared<Package>(std::bind(std::forward<Func>(func), std::forward<Args>(args)...));
        std::future<ReturnType> result = job->get_future();
        {
            std::unique_lock<std::mutex> lock(lock_);
            auto task_num = tasks_.size();
            if (stop_ || task_num > max_tasks_) throw std::runtime_error("Task queue is full/stopped!\n");

            tasks_.emplace([job]() { (*job)(); });
        }
        cv_.notify_one();
        return result;
    }
}  // namespace space::multi_thread
