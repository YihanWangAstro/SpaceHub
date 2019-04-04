
#ifndef SPACEHUB_TIMER_H
#define SPACEHUB_TIMER_H

#include <chrono>
#include <stdio.h>
#include <unistd.h>
#include <iostream>

namespace SpaceH::Tools {
    class Timer {
        using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;
    public:
        /**
         *
         */
        void start() {
            active_ = true;
            start_ = std::chrono::high_resolution_clock::now();
        }

        /**
         *
         * @return
         */
        float get_time() {
            if (active_) {
                auto now = std::chrono::high_resolution_clock::now();
                auto len = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_);
                return duration_ + static_cast<float>(len.count()) * std::chrono::milliseconds::period::num
                                   / std::chrono::milliseconds::period::den;
            } else
                return duration_;
        }

        /**
         *
         */
        void pause() {
            duration_ = get_time();
            active_ = false;
        }

        /**
         *
         */
        void reset() {
            active_ = false;
            duration_ = 0;
        }

    private:
        Time start_;
        float duration_{0};
        bool active_{false};
    };
}
#endif
