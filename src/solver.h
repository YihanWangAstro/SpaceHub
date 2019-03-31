#ifndef DYNAMICSYSTEM_H
#define DYNAMICSYSTEM_H

#include "dev-tools.h"
#include "core-computation.h"
#include <functional>

namespace SpaceH {

    template<typename ParticleSys>
    class RunArgs {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);
        using Callback = std::function<void (ParticleSys &)>;

        Scalar step_size{0};
        Scalar end_time{0};

        void pre_option(ParticleSys &partc_sys) const {
            for (auto const &opt : pre_opts_) {
                opt(partc_sys);
            }
        }

        void post_option(ParticleSys &partc_sys) const {
            for (auto const &opt : post_opts_) {
                opt(partc_sys);
            }
        }

        template<typename ...Args>
        void set_pre_option(Callback&& func, Args&& ...args) {
            pre_opts_.emplace_back(std::bind(std::forward<Callback>(func), std::placeholders::_1, std::forward<Args>(args)...));
        }

        template<typename ...Args>
        void set_post_option(Callback&& func, Args&& ...args) {
            post_opts_.emplace_back(std::bind(std::forward<Callback>(func), std::placeholders::_1, std::forward<Args>(args)...));
        }

    private:
        std::vector<Callback> pre_opts_;
        std::vector<Callback> post_opts_;
    };

    template<typename ParticSys, typename ODEiterator>
    class Solver {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);
        using RunArgs = SpaceH::RunArgs<ParticSys>;
        /* Typedef */

        template<typename STL>
        Solver(STL const &partc, Scalar t) : particles_(partc, t){}

        explicit Solver(ParticSys const& ptc, Scalar t) : particles_(ptc, t){}//more edit

        template <typename RunArgs>
        void run(const RunArgs &arg) {
            step_size_ = arg.step_size;

            if(iseq(step_size_,0.0))
                step_size_ = calc::calc_step_scale(particles_);

            Scalar end_time = arg.end_time;
            for (; particles_.time() < end_time;) {
                arg.pre_option(particles_);
                advance_one_step();
                arg.post_option(particles_);
            }
        }

        virtual ~Solver() = default; /**< @brief Default destructor, virtualize for inherent class*/
    private:
        void advance_one_step(){
            particles_.pre_iter_process();
            step_size_ = iterator_.iterate(particles_, step_size_);
            particles_.post_iter_process();
        }

        /** @brief Macro step size for ODE iterator*/
        Scalar step_size_{0.0};

        /** @brief Particle system*/
        ParticSys particles_;

        /** @brief ODE Iterator*/
        ODEiterator iterator_;
    };

}
#endif

