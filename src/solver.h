#ifndef DYNAMICSYSTEM_H
#define DYNAMICSYSTEM_H

#include "dev-tools.h"
#include "core-computation.tpp"
#include <functional>

namespace space {

    template<typename ParticleSys>
    class RunArgs {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);
        using Callback = std::function<void(ParticleSys &)>;
        using StopCall = std::function<bool(ParticleSys &)>;

        Scalar step_size{0};
        Scalar end_time{0};

        void pre_options(ParticleSys &partc_sys) const {
            for (auto const &opt : pre_opts_) {
                opt(partc_sys);
            }
        }

        void post_options(ParticleSys &partc_sys) const {
            for (auto const &opt : post_opts_) {
                opt(partc_sys);
            }
        }

        void stop_options(ParticleSys &partc_sys) const {
            for (auto const &opt : stop_opts_) {
                opt(partc_sys);
            }
        }

        bool check_stops(ParticleSys &partc_sys) const {
            for (auto const &check : stop_cond_) {
                if (check(partc_sys))
                    return true;
            }
            return false;
        }

        template<typename Func, typename ...Args>
        void add_pre_step_option(Func &&func, Args &&...args) {
            pre_opts_.emplace_back(
                    std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
        }

        template<typename Func, typename ...Args>
        void add_post_step_option(Func &&func, Args &&...args) {
            post_opts_.emplace_back(
                    std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
        }

        template<typename Func, typename ...Args>
        void add_stop_point_option(Func &&func, Args &&...args) {
            stop_opts_.emplace_back(
                    std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
        }

        template<typename Func, typename ...Args>
        void add_stop_condition(Func &&func, Args &&...args) {
            stop_cond_.emplace_back(
                    std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
        }

        template<typename Scalar>
        void add_stop_condition(Scalar end_) {
            end_time = end_;
        }

    private:
        std::vector<Callback> pre_opts_;
        std::vector<Callback> post_opts_;
        std::vector<Callback> stop_opts_;
        std::vector<StopCall> stop_cond_;
    };

    template<typename ParticSys, typename OdeIterator>
    class Solver {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);
        using RunArgs = space::RunArgs<ParticSys>;
        using Particle = typename ParticSys::Particle;

        SPACEHUB_READ_ACCESSOR(auto, particles, particles_);

        /* Typedef */
        Solver() = delete;

        Solver(Solver const &) = default;

        Solver(Solver &&) = default;

        Solver&operator=(Solver const &) = default;

        Solver&operator=(Solver &&) = default;

        template<typename STL>
        Solver(Scalar t, STL const &partc) : particles_(t, partc) {
            static_assert(is_container_v<STL>, "Only STL-like container can be used");
        }

        template<typename ...T>
        explicit Solver(Scalar t, T const &...p) : Solver(t, std::initializer_list<Particle>{p...}) {
            static_assert(calc::all(std::is_same_v<T, Particle>...), "Wrong particle type!");
        }

        void run(RunArgs const &arg) {
            step_size_ = arg.step_size;

            if (iseq(step_size_, 0.0))
                step_size_ = 0.01 * calc::calc_step_scale(particles_);

            Scalar end_time = arg.end_time;

            if (particles_.time() >= end_time)
                space::print(std::cout, "Warning: The stop time is '<=' to the start time!");

            for (; particles_.time() < end_time;) {
                if (arg.check_stops(particles_))
                    break;

                arg.pre_options(particles_);
                advance_one_step();
                arg.post_options(particles_);
            }
            arg.stop_options(particles_);
        }

        virtual ~Solver() = default; /**< @brief Default destructor, virtualize for inherent class*/
    private:
        void advance_one_step() {
            particles_.pre_iter_process();
            step_size_ = iterator_.iterate(particles_, step_size_);
            particles_.post_iter_process();
        }

        /** @brief Macro step size for ODE iterator*/
        Scalar step_size_{0.0};

        /** @brief Particle system*/
        ParticSys particles_;

        /** @brief ODE Iterator*/
        OdeIterator iterator_;
    };

}
#endif

