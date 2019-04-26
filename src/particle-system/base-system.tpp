
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "../core-computation.tpp"
#include "../dev-tools.h"
#include "../particle-system.h"
#include <type_traits>

namespace space {

    /*---------------------------------------------------------------------------*\
        Class SimpleSystem Declaration
    \*---------------------------------------------------------------------------*/
    template<typename Particles, typename Forces>
    class SimpleSystem : public ParticleSystem<SimpleSystem<Particles, Forces>> {
    public:
        //Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        //Constructors
        SimpleSystem() = delete;

        SimpleSystem(SimpleSystem const &) = default;

        SimpleSystem(SimpleSystem &&) noexcept = default;

        SimpleSystem &operator=(SimpleSystem const &) = default;

        SimpleSystem &operator=(SimpleSystem &&) noexcept = default;

        /**
         *
         * @tparam STL
         * @param t
         * @param partc
         */
        template<typename STL>
        SimpleSystem(Scalar t, STL const &partc);

        //Friend functions
        friend std::ostream &operator<<(std::ostream &os, SimpleSystem<Particles, Forces> const &ps);

        friend std::istream &operator>>(std::istream &is, SimpleSystem<Particles, Forces> &ps);

    protected:
        //CRTP implementation
        friend class ParticleSystem<SimpleSystem<Particles, Forces>>;

        /**
         *
         * @return
         */
        SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptc_.mass());

        /**
         *
         * @return
         */
        SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptc_.idn());

        /**
         *
         * @return
         */
        SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptc_.pos());

        /**
         *
         * @return
         */
        SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptc_.vel());

        /**
         *
         * @return
         */
        SPACEHUB_STD_ACCESSOR(auto, impl_time, ptc_.time());

        /**
         *
         * @return
         */
        size_t impl_number() const;

        /**
         *
         * @param dt
         */
        void impl_advance_time(Scalar dt);

        /**
         *
         * @param velocity
         * @param step_size
         */
        void impl_advance_pos(Coord const &velocity, Scalar step_size);

        /**
         *
         * @param acceleration
         * @param step_size
         */
        void impl_advance_vel(Coord const &acceleration, Scalar step_size);

        /**
         *
         * @param acceleration
         */
        void impl_evaluate_acc(Coord &acceleration) const;

        /**
         *
         * @param step_size
         */
        void impl_drift(Scalar step_size);

        /**
         *
         * @param step_size
         */
        void impl_kick(Scalar step_size);

        /**
         *
         */
        void impl_pre_iter_process();

        /**
         *
         * @tparam STL
         * @param stl
         */
        template<typename STL>
        void impl_to_linear_container(STL &stl);

        /**
         *
         * @tparam STL
         * @param stl
         */
        template<typename STL>
        void impl_load_from_linear_container(STL const &stl);

    private:
        //Private methods
        /**
         *
         */
        void eval_vel_indep_acc();

        /**
         *
         * @param step_size
         */
        void kick_pseu_vel(Scalar step_size);

        /**
         * 
         * @param step_size
         */
        void kick_real_vel(Scalar step_size);

        //Private members
        Particles ptc_;

        Forces forces_;

        Coord acc_;

        std::conditional_t<Forces::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;

        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;

        std::conditional_t<Forces::ext_vel_dep, Coord, Empty> aux_vel_;
    };

    /*---------------------------------------------------------------------------*\
        Class SimpleSystem Implementation
    \*---------------------------------------------------------------------------*/
    template<typename Particles, typename Forces>
    template<typename STL>
    SimpleSystem<Particles, Forces>::SimpleSystem(Scalar t, const STL &partc)
            : ptc_(t, partc),
              acc_(partc.size()) {
        static_assert(is_container_v<STL>, "Only STL-like container can be used");

        if constexpr (Forces::ext_vel_indep) {
            ext_vel_indep_acc_.resize(partc.size());
        }

        if constexpr (Forces::ext_vel_dep) {
            ext_vel_dep_acc_.resize(partc.size());
            aux_vel_ = ptc_.vel();
        }
    }

    template<typename Particles, typename Forces>
    template<typename STL>
    void SimpleSystem<Particles, Forces>::impl_load_from_linear_container(const STL &stl) {
        auto i = stl.begin();
        impl_time() = *i, ++i;
        load_to_coords(i, impl_pos());
        load_to_coords(i, impl_vel());
    }

    template<typename Particles, typename Forces>
    template<typename STL>
    void SimpleSystem<Particles, Forces>::impl_to_linear_container(STL &stl) {
        stl.clear();
        stl.reserve(impl_number() * 6 + 1);
        stl.emplace_back(impl_time());
        add_coords_to(stl, impl_pos());
        add_coords_to(stl, impl_vel());
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_pre_iter_process() {
        if constexpr (Forces::ext_vel_dep) {
            aux_vel_ = ptc_.vel();
        }
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_kick(Scalar step_size) {
        if constexpr (Forces::ext_vel_dep) {
            Scalar half_step = 0.5 * step_size;
            eval_vel_indep_acc();
            kick_pseu_vel(half_step);
            kick_real_vel(step_size);
            kick_pseu_vel(half_step);
        } else {
            forces_.eval_acc(ptc_, acc_);
            calc::coord_advance(ptc_.vel(), acc_, step_size);
        }
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_drift(Scalar step_size) {
        ptc_.time() += step_size;
        calc::coord_advance(ptc_.pos(), ptc_.vel(), step_size);
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_evaluate_acc(Coord &acceleration) const {
        forces_.eval_acc(ptc_, acceleration);
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_advance_vel(const Coord &acceleration, Scalar step_size) {
        calc::coord_advance(ptc_.vel(), acceleration, step_size);
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_advance_pos(const Coord &velocity, Scalar step_size) {
        calc::coord_advance(ptc_.pos(), velocity, step_size);
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::impl_advance_time(Scalar dt) {
        ptc_.time() += dt;
    }

    template<typename Particles, typename Forces>
    size_t SimpleSystem<Particles, Forces>::impl_number() const {
        return ptc_.number();
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::kick_real_vel(Scalar step_size) {
        std::swap(aux_vel_, ptc_.vel());
        forces_.eval_extra_vel_dep_acc(ptc_, ext_vel_dep_acc_);
        std::swap(aux_vel_, ptc_.vel());
        calc::coord_add(acc_, acc_, ext_vel_dep_acc_);
        calc::coord_advance(ptc_.vel(), acc_, step_size);
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::kick_pseu_vel(Scalar step_size) {
        forces_.eval_extra_vel_dep_acc(ptc_, ext_vel_dep_acc_);
        calc::coord_add(acc_, acc_, ext_vel_dep_acc_);
        calc::coord_advance(aux_vel_, acc_, step_size);
    }

    template<typename Particles, typename Forces>
    void SimpleSystem<Particles, Forces>::eval_vel_indep_acc() {
        forces_.eval_newtonian_acc(ptc_, acc_);
        if constexpr (Forces::ext_vel_dep) {
            forces_.eval_extra_vel_indep_acc(ptc_, ext_vel_indep_acc_);
            calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
        }
    }

    template<typename Particles, typename Forces>
    std::ostream &operator<<(std::ostream &os, SimpleSystem<Particles, Forces> const &ps) {
        os << ps.ptc_;
        return os;
    }

    template<typename Particles, typename Forces>
    std::istream &operator>>(std::istream &is, SimpleSystem<Particles, Forces> &ps) {
        is >> ps.ptc_;
        return is;
    }
}

#endif
