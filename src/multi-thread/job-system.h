#ifndef SPACEHUB_JOB_SYSTEM_H
#define SPACEHUB_JOB_SYSTEM_H

#include <mutex>
#include <condition_variable>
#include <tuple>
#include <fstream>
#include "safe-queue.h"

namespace SpaceH {
    namespace MultiThread {

        template<typename ...Args>
        class IOpip {
        public:
            IOpip(const char *file_name, std::ios_base::openmode mode) : file_(file_name, mode) {}

            void write_one_to_file() {
                auto data = task_.pop_front();
                file_ << data;
            }

            void read_one_from_file() {
                std::tuple<Args...> data;
                file_ >> data;
                task_.emplace_back(std::move(data));
            }

            friend IOpip& operator<<(IOpip& IO, std::tuple<Args...>&& tup){
                IO.task_.emplace_back(std::forward<decltype(tup)>(tup));
            }

            friend IOpip& operator>>(IOpip& IO, std::tuple<Args...>&& tup){
                tup = std::move(IO.task_.pop_front());
            }

        private:
            SpaceH::MultiThread::ConcurrentDeque<std::tuple<Args...>> task_;
            std::fstream file_;
        };

        class JobSystem {
        public:
        private:

        };
    }
}
#endif //SPACEHUB_JOB_SYSTEM_H
