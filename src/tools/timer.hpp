/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file timer.hpp
 *
 * Header file.
 */
#pragma once

#include <chrono>

namespace hub::tools {
    /*---------------------------------------------------------------------------*\
          Class Timer Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * High precision timer for benchmark.
     */
    class Timer {
        using Time = std::chrono::time_point<std::chrono::steady_clock>;

       public:
        /**
         * Start the timer.
         */
        void start();

        /**
         * Get the time in milli-second.
         * @return Time duration from start().
         */
        double get_time();

        /**
         * Pause the timer.
         */
        void pause();

        /**
         * Reset the timer.
         */
        void reset();

       private:
        Time start_;
        double duration_{0};
        bool active_{false};
    };
    /*---------------------------------------------------------------------------*\
          Class Timer Implementation
    \*---------------------------------------------------------------------------*/
    void Timer::start() {
        active_ = true;
        start_ = std::chrono::steady_clock::now();
    }

    double Timer::get_time() {
        if (active_) {
            auto now = std::chrono::steady_clock::now();
            auto len = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_);
            return duration_ + static_cast<double>(len.count()) * std::chrono::nanoseconds::period::num /
                                   std::chrono::nanoseconds::period::den;
        } else
            return duration_;
    }

    void Timer::pause() {
        duration_ = get_time();
        active_ = false;
    }

    void Timer::reset() {
        active_ = false;
        duration_ = 0;
    }
}  // namespace hub::tools
