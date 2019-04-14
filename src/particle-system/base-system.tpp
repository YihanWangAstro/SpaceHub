
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "../core-computation.tpp"
#include "../dev-tools.h"
#include "../particle-system.h"
#include <type_traits>

namespace space {

    template<typename Particles, typename Forces>
    class SimpleSystem : public ParticleSystem<SimpleSystem<Particles, Forces>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        SimpleSystem() = delete;

        SimpleSystem(SimpleSystem const &) = default;

        SimpleSystem(SimpleSystem &&) = default;

        SimpleSystem&operator=(SimpleSystem const &) = default;

        SimpleSystem&operator=(SimpleSystem &&) = default;

        template<typename STL>
        SimpleSystem(Scalar t, STL const &partc)
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

        friend std::ostream &operator<<(std::ostream &os, SimpleSystem const &ps) {
            os << ps.ptc_;
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SimpleSystem &ps) {
            is >> ps.ptc_;
            return is;
        }
    protected:
        friend class ParticleSystem<SimpleSystem<Particles, Forces>>;

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
            calc::coord_advance(ptc_.pos(), velocity, step_size);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar step_size) {
            calc::coord_advance(ptc_.vel(), acceleration, step_size);
        }

        void impl_evaluate_acc(Coord &acceleration) const {
            forces_.eval_acc(ptc_, acceleration);
        }

        void impl_drift(Scalar step_size) {
            ptc_.time() += step_size;
            calc::coord_advance(ptc_.pos(), ptc_.vel(), step_size);
        }

        void impl_kick(Scalar step_size) {
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

        void impl_pre_iter_process() {
            if constexpr (Forces::ext_vel_dep) {
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
            auto i = stl.begin();
            impl_time() = *i, ++i;
            load_to_coords(i, impl_pos());
            load_to_coords(i, impl_vel());
        }
    private:
        void eval_vel_indep_acc() {
            forces_.eval_newtonian_acc(ptc_, acc_);
            if constexpr (Forces::ext_vel_dep) {
                forces_.eval_extra_vel_indep_acc(ptc_, ext_vel_indep_acc_);
                calc::coord_add(acc_, acc_, ext_vel_indep_acc_);
            }
        }

        void kick_pseu_vel(Scalar step_size) {
            forces_.eval_extra_vel_dep_acc(ptc_, ext_vel_dep_acc_);
            calc::coord_add(acc_, acc_, ext_vel_dep_acc_);
            calc::coord_advance(aux_vel_, acc_, step_size);
        }

        void kick_real_vel(Scalar step_size) {
            std::swap(aux_vel_, ptc_.vel());
            forces_.eval_extra_vel_dep_acc(ptc_, ext_vel_dep_acc_);
            std::swap(aux_vel_, ptc_.vel());
            calc::coord_add(acc_, acc_, ext_vel_dep_acc_);
            calc::coord_advance(ptc_.vel(), acc_, step_size);
        }

        Particles ptc_;
        Forces forces_;
        Coord acc_;

        std::conditional_t <Forces::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;
        std::conditional_t <Forces::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;
        std::conditional_t <Forces::ext_vel_dep, Coord, Empty> aux_vel_;
    };
}

#endif
