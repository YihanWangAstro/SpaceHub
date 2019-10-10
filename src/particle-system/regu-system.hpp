
#ifndef REGUPARTICLESYSTEM_H
#define REGUPARTICLESYSTEM_H

#include <type_traits>
#include "../core-computation.hpp"
#include "../dev-tools.hpp"
#include "../particle-system.hpp"
#include "../accelerations.hpp"

namespace space {
  enum class ReguType {
    logH, TTL, none
  };

/*---------------------------------------------------------------------------*\
    Class Regularization Declaration
\*---------------------------------------------------------------------------*/
  template<typename Scalar, ReguType Type = ReguType::logH>
  class Regularization {
  public:
    // Constructors
    template<typename Particles>
    explicit Regularization(Particles const &partc);

    // Public methods
    SPACEHUB_STD_ACCESSOR(auto, omega, omega_);

    SPACEHUB_STD_ACCESSOR(auto, bindE, bindE_);

    template<typename Particles>
    Scalar eval_pos_phy_time(Particles const &partc, Scalar step_size) const;

    template<typename Particles>
    Scalar eval_vel_phy_time(Particles const &partc, Scalar step_size) const;

  private:
    // Private methods
    template<typename Particles>
    inline Scalar capital_omega(Particles const &partc) const;

    // Private members
    Scalar omega_;

    Scalar bindE_;
  };

/*---------------------------------------------------------------------------*\
    Class RegularizedSystem Declaration
\*---------------------------------------------------------------------------*/

/**
 * @brief Regularized particle System.
 *
 * Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
 * http://iopscience.iop.org/article/10.1086/301102/meta and
 * https://link.springer.com/article/10.1023%2FA%3A1021149313347.
 * @tparam Particles
 * @tparam Interactions
 */
  template<typename Particles, typename Interactions, ReguType RegType>
  class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Interactions, RegType>> {
  public:
    // Type members
    SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

    using Base = ParticleSystem<RegularizedSystem<Particles, Interactions, RegType>>;

    using Particle = typename Particles::Particle;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(RegularizedSystem, delete, default, default, default, default);

    template<typename STL>
    RegularizedSystem(Scalar t, STL const &particle_set);

    // Static members
    static constexpr ReguType regu_type{RegType};

    // Public methods
    SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

    SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

    // Friend functions
    template<typename P, typename F, ReguType R>
    friend std::ostream &operator<<(std::ostream &os, RegularizedSystem<P, F, R> const &ps);

    template<typename P, typename F, ReguType R>
    friend std::istream &operator>>(std::istream &is, RegularizedSystem<P, F, R> &ps);

    CRTP_impl :
    // CRTP implementation
    SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptcl_.mass());

    SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptcl_.idn());

    SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptcl_.pos());

    SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptcl_.vel());

    SPACEHUB_STD_ACCESSOR(auto, impl_time, ptcl_.time());

    size_t impl_number() const;

    void impl_advance_time(Scalar step_size);

    void impl_advance_pos(Coord const &velocity, Scalar step_size);

    void impl_advance_vel(Coord const &acceleration, Scalar step_size);

    void impl_evaluate_acc(Coord &acceleration) const;

    void impl_drift(Scalar step_size);

    void impl_kick(Scalar step_size);

    void impl_pre_iter_process();

    template<typename STL>
    void impl_to_linear_container(STL &stl);

    template<typename STL>
    void impl_load_from_linear_container(STL const &stl);

  private:
    // Private methods
    void eval_vel_indep_acc();

    void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar phy_time);

    void advance_bindE(Coord const &velocity, Coord const &d_bindE_dr, Scalar phy_time);

    void kick_pseu_vel(Scalar phy_time);

    void kick_real_vel(Scalar phy_time);

    // Private members
    Particles ptcl_;
    Interactions interactions_;
    Accelerations<Interactions, Coord> accels_;
    Regularization<Scalar, RegType> regu_;
    std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> aux_vel_;
  };

/*---------------------------------------------------------------------------*\
    Class RegularizedSystem Implementation
\*---------------------------------------------------------------------------*/
  template<typename Particles, typename Interactions, ReguType RegType>
  template<typename STL>
  RegularizedSystem<Particles, Interactions, RegType>::RegularizedSystem(Scalar t, const STL &particle_set)
          : ptcl_(t, particle_set), accels_(particle_set.size()), regu_(ptcl_) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
    if constexpr (Interactions::ext_vel_dep) {
      aux_vel_ = ptcl_.vel();
    }
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  std::istream &operator>>(std::istream &is, RegularizedSystem<Particles, Interactions, RegType> &ps) {
    is >> ps.ptcl_;
    return is;
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  std::ostream &operator<<(std::ostream &os, const RegularizedSystem<Particles, Interactions, RegType> &ps) {
    os << ps.ptcl_;
    return os;
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  size_t RegularizedSystem<Particles, Interactions, RegType>::impl_number() const {
    return ptcl_.number();
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_advance_time(Scalar step_size) {
    Scalar phy_time = regu_.eval_pos_phy_time(ptcl_, step_size);
    ptcl_.time() += phy_time;
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_advance_pos(const Coord &velocity, Scalar step_size) {
    Scalar phy_time = regu_.eval_pos_phy_time(ptcl_, step_size);
    calc::coord_advance(ptcl_.pos(), velocity, phy_time);
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void
  RegularizedSystem<Particles, Interactions, RegType>::impl_advance_vel(const Coord &acceleration, Scalar step_size) {
    Scalar phy_time = regu_.eval_vel_phy_time(ptcl_, step_size);
    calc::coord_advance(ptcl_.vel(), acceleration, phy_time);
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_evaluate_acc(Coord &acceleration) const {
    interactions_.eval_acc(ptcl_, acceleration);
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_drift(Scalar step_size) {
    Scalar phy_time = regu_.eval_pos_phy_time(ptcl_, step_size);
    calc::coord_advance(ptcl_.pos(), ptcl_.vel(), phy_time);
    ptcl_.time() += phy_time;
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_kick(Scalar step_size) {
    Scalar phy_time = regu_.eval_vel_phy_time(ptcl_, step_size);
    Scalar half_time = 0.5 * phy_time;

    eval_vel_indep_acc();

    if constexpr (Interactions::ext_vel_dep) {
      kick_pseu_vel(half_time);
      kick_real_vel(phy_time);
      kick_pseu_vel(half_time);
    } else {
      calc::coord_advance(ptcl_.vel(), accels_.tot_vel_indep_acc(), half_time);
      advance_omega(ptcl_.vel(), accels_.newtonian_acc(), phy_time);
      if constexpr (Interactions::ext_vel_indep) {
        advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), phy_time);
      }
      calc::coord_advance(ptcl_.vel(), accels_.tot_vel_indep_acc(), half_time);
      /*advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
      if constexpr (Interactions::ext_vel_indep) {
        advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
      }
      calc::coord_advance(ptcl_.vel(), accels_.tot_vel_indep_acc(), phy_time);
      advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
      if constexpr (Interactions::ext_vel_indep) {
        advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
      }*/
    }
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_pre_iter_process() {
    if constexpr (Interactions::ext_vel_dep) {
      aux_vel_ = ptcl_.vel();
    }
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  template<typename STL>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_to_linear_container(STL &stl) {
    stl.clear();
    stl.reserve(impl_number() * 6 + 3);
    stl.emplace_back(impl_time());
    stl.emplace_back(omega());
    stl.emplace_back(bindE());
    add_coords_to(stl, impl_pos());
    add_coords_to(stl, impl_vel());
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  template<typename STL>
  void RegularizedSystem<Particles, Interactions, RegType>::impl_load_from_linear_container(const STL &stl) {
    auto begin = stl.begin();
    impl_time() = *begin;
    omega() = *(begin + 1);
    bindE() = *(begin + 2);

    size_t len = impl_number() * 3;

    auto pos_begin = begin + 3;
    auto pos_end = pos_begin + len;
    auto vel_begin = pos_end;
    auto vel_end = vel_begin + len;
    load_to_coords(pos_begin, pos_end, impl_pos());
    load_to_coords(vel_begin, vel_end, impl_vel());
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::eval_vel_indep_acc() {
    interactions_.eval_newtonian_acc(ptcl_, accels_.newtonian_acc());
    accels_.tot_vel_indep_acc() = accels_.newtonian_acc();
    if constexpr (Interactions::ext_vel_indep) {
      interactions_.eval_extra_vel_indep_acc(ptcl_, accels_.ext_vel_indep_acc());
      calc::coord_add(accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc(), accels_.newtonain_acc());
    }
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void
  RegularizedSystem<Particles, Interactions, RegType>::advance_omega(const Coord &velocity, const Coord &d_omega_dr,
                                                                     Scalar phy_time) {
    Scalar d_omega = calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_omega_dr);
    regu_.omega() += d_omega * phy_time;
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void
  RegularizedSystem<Particles, Interactions, RegType>::advance_bindE(const Coord &velocity, const Coord &d_bindE_dr,
                                                                     Scalar phy_time) {
    Scalar d_bindE = -calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_bindE_dr);
    regu_.bindE() += d_bindE * phy_time;
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
    interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
    calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
    calc::coord_advance(aux_vel_, accels_.acc(), phy_time);
  }

  template<typename Particles, typename Interactions, ReguType RegType>
  void RegularizedSystem<Particles, Interactions, RegType>::kick_real_vel(Scalar phy_time) {
    std::swap(aux_vel_, ptcl_.vel());
    interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
    std::swap(aux_vel_, ptcl_.vel());

    calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
    calc::coord_advance(ptcl_.vel(), accels_.acc(), phy_time);

    advance_omega(aux_vel_, accels_.newtonian_acc(), phy_time);

    if constexpr (Interactions::ext_vel_indep) {
      calc::coord_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
      advance_bindE(aux_vel_, accels_.acc(), phy_time);
    } else {
      advance_bindE(aux_vel_, accels_.ext_vel_dep_acc(), phy_time);
    }
  }

/*---------------------------------------------------------------------------*\
    Class Regularization Implementation
\*---------------------------------------------------------------------------*/
  template<typename Scalar, ReguType Type>
  template<typename Particles>
  Regularization<Scalar, Type>::Regularization(const Particles &partc) {
    omega_ = capital_omega(partc);
    bindE_ = -calc::calc_total_energy(partc);
  }

  template<typename Scalar, ReguType Type>
  template<typename Particles>
  Scalar Regularization<Scalar, Type>::eval_pos_phy_time(const Particles &partc, Scalar step_size) const {
    if constexpr (Type == ReguType::logH) {
      return step_size / (bindE_ + calc::calc_kinetic_energy(partc));
    } else if constexpr (Type == ReguType::TTL) {
      return step_size / omega_;
    } else if constexpr (Type == ReguType::none) {
      return step_size;
    } else {
      spacehub_abort("Undefined regularization type!");
    }
  }

  template<typename Scalar, ReguType Type>
  template<typename Particles>
  Scalar Regularization<Scalar, Type>::eval_vel_phy_time(const Particles &partc, Scalar step_size) const {
    if constexpr (Type == ReguType::logH) {
      return step_size / -calc::calc_potential_energy(partc);
    } else if constexpr (Type == ReguType::TTL) {
      return step_size / capital_omega(partc);
    } else if constexpr (Type == ReguType::none) {
      return step_size;
    } else {
      spacehub_abort("Undefined regularization type!");
    }
  }

  template<typename Scalar, ReguType Type>
  template<typename Particles>
  Scalar Regularization<Scalar, Type>::capital_omega(const Particles &partc) const {
    return -calc::calc_potential_energy(partc);
  }
}  // namespace space

#endif
