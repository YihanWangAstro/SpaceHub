
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "../core-computation.tpp"
#include "../dev-tools.h"
#include "../particle-system.h"


namespace SpaceH {

    template<typename Particles, typename Interactions>
    class SimpleSystem : public ParticleSystem<SimpleSystem<Particles, Interactions>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        SimpleSystem() = delete;

        SimpleSystem(SimpleSystem const &) = default;

        SimpleSystem(SimpleSystem &&) = default;

        SimpleSystem&operator=(SimpleSystem const &) = default;

        SimpleSystem&operator=(SimpleSystem &&) = default;

        template<typename STL>
        SimpleSystem(Scalar t, STL const &partc) : ptc_(t, partc), acc_(partc.size()) {
            if constexpr (Interactions::has_extra_vel_indep_acc) {
                extra_vel_indep_acc_.resize(partc.size());
            }

            if constexpr (Interactions::has_extra_vel_dep_acc) {
                extra_vel_dep_acc_.resize(partc.size());
                aux_vel_ = ptc_.vel();
            }
        }

        friend std::ostream &operator<<(std::ostream &os, SimpleSystem const &ps) {
            os << ps.ptc_;
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SimpleSystem &ps) {
            is >> ps.ptc_;
            return is;
        }

    protected:
        friend class ParticleSystem<SimpleSystem<Particles, Interactions>>;

        SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptc_.mass());

        SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptc_.idn());

        SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptc_.pos());

        SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptc_.vel());

        SPACEHUB_STD_ACCESSOR(auto, impl_time, ptc_.time());

        size_t impl_number() const {
            return ptc_.number();
        }

        void impl_advance_time(Scalar dt) {
            ptc_.time() += dt;
        }

        void impl_advance_pos(Coord const &velocity, Scalar step_size) {
            Calc::coord_advance(ptc_.pos(), velocity, step_size);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar step_size) {
            Calc::coord_advance(ptc_.vel(), acceleration, step_size);
        }

        void impl_evaluate_acc(Coord &acceleration) const {
            eom_.eval_acc(ptc_, acceleration);
        }

        void impl_drift(Scalar step_size) {
            ptc_.time() += step_size;
            Calc::coord_advance(ptc_.pos(), ptc_.vel(), step_size);
        }

        void impl_kick(Scalar step_size) {
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                Scalar half_step = 0.5 * step_size;
                eval_vel_indep_acc();
                kick_pseu_vel(half_step);
                kick_real_vel(step_size);
                kick_pseu_vel(half_step);
            } else {
                eom_.eval_acc(ptc_, acc_);
                Calc::coord_advance(ptc_.vel(), acc_, step_size);
            }
        }

        void impl_pre_iter_process() {
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                aux_vel_ = ptc_.vel();
            }
        }

        template <typename STL>
        void impl_to_linear_container(STL& stl){
            stl.clear();
            stl.reserve(impl_number()*6 +1);
            stl.emplace_back(impl_time());
            add_coords_to(stl, impl_pos());
            add_coords_to(stl, impl_vel());
        }

        template <typename STL>
        void impl_load_from_linear_container(STL const& stl){
            size_t i = 0;
            impl_time() = stl[i++];
            load_to_coords(stl, i, impl_pos());
            load_to_coords(stl, i, impl_vel());
        }
    private:
        void eval_vel_indep_acc() {
            eom_.eval_newtonian_acc(ptc_, acc_);
            if constexpr (Interactions::has_extra_vel_dep_acc) {
                eom_.eval_extra_vel_indep_acc(ptc_, extra_vel_indep_acc_);
                Calc::coord_add(acc_, acc_, extra_vel_indep_acc_);
            }
        }

        void kick_pseu_vel(Scalar step_size) {
            eom_.eval_extra_vel_dep_acc(ptc_, extra_vel_dep_acc_);
            Calc::coord_add(acc_, acc_, extra_vel_dep_acc_);
            Calc::coord_advance(aux_vel_, acc_, step_size);
        }

        void kick_real_vel(Scalar step_size) {
            std::swap(aux_vel_, ptc_.vel());
            eom_.eval_extra_vel_dep_acc(ptc_, extra_vel_dep_acc_);
            std::swap(aux_vel_, ptc_.vel());
            Calc::coord_add(acc_, acc_, extra_vel_dep_acc_);
            Calc::coord_advance(ptc_.vel(), acc_, step_size);
        }

        Particles ptc_;
        Interactions eom_;
        Coord acc_{0};

        Coord extra_vel_indep_acc_{0};
        Coord extra_vel_dep_acc_{0};
        Coord aux_vel_{0};
    };
}

#endif
