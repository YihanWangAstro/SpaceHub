
#ifndef SPACEHUB_BASE_SYSTEM_HPP
#define SPACEHUB_BASE_SYSTEM_HPP

#include <type_traits>
#include "../core-computation.hpp"
#include "../particle-system.hpp"
#include "../accelerations.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class SimpleSystem Declaration
\*---------------------------------------------------------------------------*/
  template<typename Particles, typename Interactions>
  class SimpleSystem : public ParticleSystem<SimpleSystem<Particles, Interactions>> {
  public:
    // Type members
    SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

    using Base = ParticleSystem<SimpleSystem<Particles, Interactions>>;

    using Particle = typename Particles::Particle;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(SimpleSystem, delete, default, default, default, default);

    /**
     *
     * @tparam STL
     * @param time
     * @param particle_set
     */
    template<typename STL>
    SimpleSystem(Scalar time, STL const &particle_set);

    // Friend functions
    template<typename P, typename F>
    friend std::ostream &operator<<(std::ostream &os, SimpleSystem<P, F> const &ps);

    template<typename P, typename F>
    friend std::istream &operator>>(std::istream &is, SimpleSystem<P, F> &ps);

    CRTP_IMPL :
    // CRTP implementation

    /**
     *
     * @return
     */
    SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptcl_.mass());

    /**
     *
     * @return
     */
    SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptcl_.idn());

    /**
     *
     * @return
     */
    SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptcl_.pos());

    /**
     *
     * @return
     */
    SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptcl_.vel());

    /**
     *
     * @return
     */
    SPACEHUB_STD_ACCESSOR(auto, impl_time, ptcl_.time());

    /**
     *
     * @return
     */
    [[nodiscard]] size_t impl_number() const;

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
     * @param stl_ranges
     */
    template<typename STL>
    void impl_to_linear_container(STL &stl_ranges);

    /**
     *
     * @tparam STL
     * @param stl_ranges
     */
    template<typename STL>
    void impl_load_from_linear_container(STL const &stl_ranges);

  private:
    // Private methods
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

    // Private members
    Particles ptcl_;

    Interactions interactions_;

    Accelerations <Interactions, Coord> accels_;

    std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> aux_vel_;
  };

/*---------------------------------------------------------------------------*\
    Class SimpleSystem Implementation
\*---------------------------------------------------------------------------*/
  template<typename Particles, typename Interactions>
  template<typename STL>
  SimpleSystem<Particles, Interactions>::SimpleSystem(Scalar time, const STL &particle_set) : ptcl_(time, particle_set),
                                                                                              accels_(particle_set.size()) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
    if constexpr (Interactions::ext_vel_dep) {
      aux_vel_ = ptcl_.vel();
    }
  }

  template<typename Particles, typename Interactions>
  template<typename STL>
  void SimpleSystem<Particles, Interactions>::impl_load_from_linear_container(const STL &stl_ranges) {
    auto begin = stl_ranges.begin();
    impl_time() = *begin;
    size_t len = impl_number() * 3;
    auto pos_begin = begin + 1;
    auto pos_end = pos_begin + len;
    auto vel_begin = pos_end;
    auto vel_end = vel_begin + len;
    load_to_coords(pos_begin, pos_end, impl_pos());
    load_to_coords(vel_begin, vel_end, impl_vel());
  }

  template<typename Particles, typename Interactions>
  template<typename STL>
  void SimpleSystem<Particles, Interactions>::impl_to_linear_container(STL &stl_ranges) {
    stl_ranges.clear();
    stl_ranges.reserve(impl_number() * 6 + 1);
    stl_ranges.emplace_back(impl_time());
    add_coords_to(stl_ranges, impl_pos());
    add_coords_to(stl_ranges, impl_vel());
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_pre_iter_process() {
    if constexpr (Interactions::ext_vel_dep) {
      aux_vel_ = ptcl_.vel();
    }
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_kick(Scalar step_size) {
    if constexpr (Interactions::ext_vel_dep) {
      Scalar half_step = 0.5 * step_size;
      eval_vel_indep_acc();
      kick_pseu_vel(half_step);
      kick_real_vel(step_size);
      kick_pseu_vel(half_step);
    } else {
      interactions_.eval_acc(ptcl_, accels_.acc());
      calc::coord_advance(ptcl_.vel(), accels_.acc(), step_size);
    }
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_drift(Scalar step_size) {
    ptcl_.time() += step_size;
    calc::coord_advance(ptcl_.pos(), ptcl_.vel(), step_size);
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_evaluate_acc(Coord &acceleration) const {
    interactions_.eval_acc(ptcl_, acceleration);
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_advance_vel(const Coord &acceleration, Scalar step_size) {
    calc::coord_advance(ptcl_.vel(), acceleration, step_size);
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_advance_pos(const Coord &velocity, Scalar step_size) {
    calc::coord_advance(ptcl_.pos(), velocity, step_size);
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::impl_advance_time(Scalar dt) {
    ptcl_.time() += dt;
  }

  template<typename Particles, typename Interactions>
  size_t SimpleSystem<Particles, Interactions>::impl_number() const {
    return ptcl_.number();
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::kick_real_vel(Scalar step_size) {
    std::swap(aux_vel_, ptcl_.vel());
    interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
    std::swap(aux_vel_, ptcl_.vel());
    calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
    calc::coord_advance(ptcl_.vel(), accels_.acc(), step_size);
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::kick_pseu_vel(Scalar step_size) {
    interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
    calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
    calc::coord_advance(aux_vel_, accels_.acc(), step_size);
  }

  template<typename Particles, typename Interactions>
  void SimpleSystem<Particles, Interactions>::eval_vel_indep_acc() {
    interactions_.eval_newtonian_acc(ptcl_, accels_.tot_vel_indep_acc());
    if constexpr (Interactions::ext_vel_indep) {
      interactions_.eval_extra_vel_indep_acc(ptcl_, accels_.ext_vel_indep_acc());
      calc::coord_add(accels_.tot_vel_indep_acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc());
    }
  }

  template<typename Particles, typename Interactions>
  std::ostream &operator<<(std::ostream &os, SimpleSystem<Particles, Interactions> const &ps) {
    os << ps.ptcl_;
    return os;
  }

  template<typename Particles, typename Interactions>
  std::istream &operator>>(std::istream &is, SimpleSystem<Particles, Interactions> &ps) {
    is >> ps.ptcl_;
    return is;
  }
}  // namespace space

#endif
