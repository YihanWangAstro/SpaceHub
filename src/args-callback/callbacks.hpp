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
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 *
 * @file callbacks.hpp
 *
 * Header file.
 */
#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>

#include "../dev-tools.hpp"

/**
 * @namespace space::callback
 * Documentation for callback
 */
namespace space::callback {

    /*---------------------------------------------------------------------------*\
     Class TimeSlice Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Create a wrapper on callable object(function pointer, functor, lambda) of
     * which accepts one parameter. The type of the parameter must have method
     * `time()`. The wrapped callable object will be invoked for equal spacing time
     * interval for time provided by time(). This class is basically used to
     * implement equal time operation in simulations. For example, one can provide a
     * printer `[](auto&p){std::cout << p << std::endl;}` as the pre_step_operation
     * in RunArgs to output the state of the integrated system. This printer then
     * will be invoked before every step integration. The output might be very dense
     * sometimes, thus outputting the result of every step is somewhat heavy. If one
     * wraps the printer with TimeSlice, `TimeSlice([](auto&p){std::cout << p <<
     * '\n'}, 0, 100, 10)`, then the output will be performed only at
     * p.time()=[0,100/10, 2*100/10, 3*100/10,...100].
     * @tparam Operation Callable object.
     */
    template <typename Operation, typename Scalar>
    class TimeSlice {
       public:
        SPACEHUB_MAKE_CONSTRUCTORS(TimeSlice, default, default, default, default, default);

        /**
         * Constructor of the time slice
         * @param[in] opt Callable object
         * @param[in] start The start time of the time slice
         * @param[in] end The end time of the time slice
         * @param[in] opt_num The bin number of the time slice. i.e time interval =
         * (end-start)/opt_num
         */
        TimeSlice(Operation const& opt, Scalar start, Scalar end, size_t opt_num = 5000);

        /**
         * Callable interface.
         * @tparam ParticleSys Any type provides method `time()`
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         */
        template <typename ParticleSys>
        inline auto operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
            std::is_same_v<void, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, void>;

        /**
         * Callable interface.
         * @tparam ParticleSys Any type provides method `time()`
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         * @return auto bool
         */
        template <typename ParticleSys>
        inline auto operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
            std::is_same_v<bool, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, bool>;

        /**
         * Reset the slice parameters.
         * @param[in] start The start time of the time slice.
         * @param[in] end The end time of the time slice.
         * @param[in] opt_num The bin number of the time slice. i.e time interval =
         * (end-start)/opt_num
         */
        void reset_slice_params(Scalar start, Scalar end, size_t opt_num = 5000);

        Operation operation() { return opt_; };

       private:
        Operation opt_;
        Scalar opt_time_{0};
        Scalar end_time_{0};
        Scalar opt_interval_{0};
    };

    /*---------------------------------------------------------------------------*\
     Class StepSlice Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Create a wrapper on callable object(function pointer, functor, lambda) of
     * which accepts one parameter. The wrapped callable object will be invoked for
     * equal spacing steps. This class is basically used to implement equal step
     * operation in simulations. For example, one can provide a printer
     * `[](auto&p){std::cout << p << std::endl;}` as the pre_step_operation in
     * RunArgs to output the state of the integrated system. This printer then will
     * be invoked before every step integration. The output might be very dense
     * sometimes, thus outputting the result of every step is somewhat heavy. If one
     * wraps the printer with StepSlice, `StepSlice([](auto&p){std::cout << p <<
     * '\n'}, 10)`, then the output will be performed every 10 steps.
     * @tparam Operation Callable object.
     */
    template <typename Operation>
    class StepSlice {
       public:
        SPACEHUB_MAKE_CONSTRUCTORS(StepSlice, default, default, default, default, default);

        /**
         * Constructor of step slice.
         * @param[in] opt Callable object.
         * @param[in] step_interval Step interval.
         */
        explicit StepSlice(Operation const& opt, size_t step_interval = 1);

        /**
         * Callable interface.
         * @tparam ParticleSys Any type used as call back parameter.
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         */
        template <typename ParticleSys>
        inline auto operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
            std::is_same_v<void, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, void>;

        /**
         * Callable interface.
         * @tparam ParticleSys Any type used as call back parameter.
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         * @return auto bool.
         */
        template <typename ParticleSys>
        inline auto operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
            std::is_same_v<bool, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, bool>;

        /**
         * Reset the slice parameters.
         * @param[in] step_interval Step interval.
         */
        void reset_slice_params(size_t step_interval);

        Operation operation() { return opt_; };

       private:
        Operation opt_;
        size_t step_{0};
        size_t step_interval_{0};
    };

    /*---------------------------------------------------------------------------*\
        Class LogTimeSlice Declaration
       \*---------------------------------------------------------------------------*/
    /**

     * @tparam Operation Callable object.
     */
    template <typename Operation, typename Scalar>
    class LogTimeSlice {
       public:
        SPACEHUB_MAKE_CONSTRUCTORS(LogTimeSlice, default, default, default, default, default);

        /**
         * Constructor of the time slice
         * @param[in] opt Callable object
         * @param[in] start The start time of the time slice
         * @param[in] end The end time of the time slice
         * @param[in] opt_num The bin number of the time slice. i.e time interval =
         * (end-start)/opt_num
         */
        LogTimeSlice(Operation const& opt, Scalar start, Scalar end, size_t opt_num = 5000);

        /**
         * Callable interface.
         * @tparam ParticleSys Any type provides method `time()`
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         */
        template <typename ParticleSys>
        inline auto operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
            std::is_same_v<void, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, void>;

        /**
         * Callable interface.
         * @tparam ParticleSys Any type provides method `time()`
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         * @return auto bool
         */
        template <typename ParticleSys>
        inline auto operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
            std::is_same_v<bool, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, bool>;

        /**
         * Reset the slice parameters.
         * @param[in] start The start time of the time slice.
         * @param[in] end The end time of the time slice.
         * @param[in] opt_num The bin number of the time slice. i.e time interval =
         * (end-start)/opt_num
         */
        void reset_slice_params(Scalar start, Scalar end, size_t opt_num = 5000);

        Operation operation() { return opt_; };

       private:
        Operation opt_;
        Scalar opt_time_{0};
        Scalar start_time_{0};
        Scalar end_time_{0};
        Scalar opt_interval_{1};
    };
    /*---------------------------------------------------------------------------*\
     Class DefaultWriter Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Default output writer for RunArgs. This class serves as a callable callback object
     * to output data to a file stream.
     */
    class DefaultWriter {
       public:
        /**
         * Constructor of the output writer.
         * @param file_name The file name of the output file stream.
         */
        explicit DefaultWriter(std::string const& file_name);

        /**
         * Callable interface.
         * @tparam ParticleSys Any type provides method `time()`
         * @param[in,out] ptc Input parameter.
         * @param[in] step_size The step size of the integration.
         */
        template <typename ParticleSys>
        inline void operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size);

        template <typename T>
        friend DefaultWriter& operator<<(DefaultWriter& wtr, T const& d);

       private:
        std::shared_ptr<std::ofstream> fstream_;
    };

    /*---------------------------------------------------------------------------*\
     Class TimeSlice Definition
    \*---------------------------------------------------------------------------*/
    template <typename Operation, typename Scalar>
    TimeSlice<Operation, Scalar>::TimeSlice(const Operation& opt, Scalar start, Scalar end, size_t opt_num)
        : opt_{opt}, opt_time_{start}, end_time_{end}, opt_interval_{(end - start) / opt_num} {
        assert(end >= start);
    }

    template <typename Operation, typename Scalar>
    template <typename ParticleSys>
    auto TimeSlice<Operation, Scalar>::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size)
        -> std::enable_if_t<
            std::is_same_v<void, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, void> {
        auto t = ptc.time();
        if (t >= static_cast<Scalar>(opt_time_) && t <= end_time_) {
            opt_time_ += opt_interval_;
            opt_(ptc, step_size);
        }
    }

    template <typename Operation, typename Scalar>
    template <typename ParticleSys>
    auto TimeSlice<Operation, Scalar>::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size)
        -> std::enable_if_t<
            std::is_same_v<bool, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, bool> {
        auto t = ptc.time();
        if (ptc.time() >= static_cast<Scalar>(opt_time_) && t <= end_time_) {
            opt_time_ += opt_interval_;
            return opt_(ptc, step_size);
        } else {
            return false;
        }
    }

    template <typename Operation, typename Scalar>
    void TimeSlice<Operation, Scalar>::reset_slice_params(Scalar start, Scalar end, size_t opt_num) {
        opt_time_ = start;
        end_time_ = end;
        opt_interval_ = (end - start) / opt_num;
    }

    /*---------------------------------------------------------------------------*\
      Class StepSlice Definition
    \*---------------------------------------------------------------------------*/
    template <typename Operation>
    StepSlice<Operation>::StepSlice(const Operation& opt, size_t step_interval)
        : opt_{opt}, step_interval_{step_interval} {}

    template <typename Operation>
    template <typename ParticleSys>
    auto StepSlice<Operation>::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
        std::is_same_v<void, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, void> {
        if (step_ % step_interval_ == 0) {
            opt_(ptc, step_size);
        }
        step_++;
    }

    template <typename Operation>
    template <typename ParticleSys>
    auto StepSlice<Operation>::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) -> std::enable_if_t<
        std::is_same_v<bool, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, bool> {
        if (step_ % step_interval_ == 0) {
            step_++;
            return opt_(ptc, step_size);
        } else {
            step_++;
            return false;
        }
    }

    template <typename Operation>
    void StepSlice<Operation>::reset_slice_params(size_t step_interval) {
        step_ = 0;
        step_interval_ = step_interval;
    }

    /*---------------------------------------------------------------------------*\
        Class LogTimeSlice Definition //TODO
    \*---------------------------------------------------------------------------*/
    template <typename Operation, typename Scalar>
    LogTimeSlice<Operation, Scalar>::LogTimeSlice(const Operation& opt, Scalar start, Scalar end, size_t opt_num)
        : opt_{opt}, opt_time_{start}, start_time_{start}, end_time_{end}, opt_interval_{log(end - start) / opt_num} {
        assert(end >= start);
    }

    template <typename Operation, typename Scalar>
    template <typename ParticleSys>
    auto LogTimeSlice<Operation, Scalar>::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size)
        -> std::enable_if_t<
            std::is_same_v<void, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, void> {
        auto t = ptc.time();
        if (t >= static_cast<Scalar>(opt_time_) && t <= end_time_) {
            opt_time_ = (opt_time_ - start_time_) * opt_interval_ + start_time_;
            opt_(ptc, step_size);
        }
    }

    template <typename Operation, typename Scalar>
    template <typename ParticleSys>
    auto LogTimeSlice<Operation, Scalar>::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size)
        -> std::enable_if_t<
            std::is_same_v<bool, std::result_of_t<Operation(ParticleSys&, typename ParticleSys::Scalar)>>, bool> {

        auto t = ptc.time();
        if (ptc.time() >= static_cast<Scalar>(opt_time_) && t <= end_time_) {
            opt_time_ += opt_interval_;
            return opt_(ptc, step_size);
        } else {
            return false;
        }
    }

    template <typename Operation, typename Scalar>
    void LogTimeSlice<Operation, Scalar>::reset_slice_params(Scalar start, Scalar end, size_t opt_num) {
        opt_time_ = start;
        end_time_ = end;
        opt_interval_ = (end - start) / opt_num;
    }

    /*---------------------------------------------------------------------------*\
       Class DefaultWriter Definition
    \*---------------------------------------------------------------------------*/
    DefaultWriter::DefaultWriter(std::string const& file_name) : fstream_{std::make_shared<std::ofstream>(file_name)} {
        if (!fstream_->is_open()) {
            spacehub_abort("Fail to open the file " + file_name);
        } else {
            (*fstream_) << std::scientific << std::setprecision(16);
        }
    }

    template <typename ParticleSys>
    void DefaultWriter::operator()(ParticleSys& ptc, typename ParticleSys::Scalar step_size) {
        if (ptc.time() == 0) [[unlikely]] {
            *fstream_ << ptc.column_names() << '\n';
        }
        *fstream_ << ptc << std::endl;  //'\n';
    }

    template <typename T>
    DefaultWriter& operator<<(DefaultWriter& wtr, const T& d) {
        (*wtr.fstream_) << d;
        return wtr;
    }
}  // namespace space::callback
