//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_NEWTONIAN_HPP
#define SPACEHUB_NEWTONIAN_HPP

#include "interaction.hpp"
/**
 * @namespace space::interactions
 * Documentation for space
 */
namespace space::interactions {
 /*---------------------------------------------------------------------------*\
      Class NewtonianGrav Declaration
\*---------------------------------------------------------------------------*/
  /**
   *
   */
  class NewtonianGrav : public Interactions<NewtonianGrav> {
  public:
    //Type members
    using Base = Interactions<NewtonianGrav>;
    CRTP_IMPL:

    //CRTP implementation
    template<typename Particles>
    void impl_eval_acc(Particles const &particles, typename Particles::Coord &acceleration) const;

    template<typename Particles>
    void impl_eval_newtonian_acc(Particles const &particles, typename Particles::Coord &acceleration) const;

  private:
    CREATE_METHOD_CHECK(chain_pos);

    CREATE_METHOD_CHECK(index);
  };

/*---------------------------------------------------------------------------*\
      Class NewtonianGrav Implememtation
\*---------------------------------------------------------------------------*/
  template<typename Particles>
  void NewtonianGrav::impl_eval_acc(const Particles &particles, typename Particles::Coord &acceleration) const {
    impl_eval_newtonian_acc(particles, acceleration);
  }

  template<typename Particles>
  void NewtonianGrav::impl_eval_newtonian_acc(const Particles &particles, typename Particles::Coord &acceleration) const {
    size_t num = particles.number();
    auto &px = particles.pos().x;
    auto &py = particles.pos().y;
    auto &pz = particles.pos().z;
    auto &m = particles.mass();

    calc::set_arrays_zero(acceleration.x, acceleration.y, acceleration.z);

    auto force = [&](auto dx, auto dy, auto dz, auto i, auto j) {
      auto r = sqrt(dx * dx + dy * dy + dz * dz);
      auto rr3 = 1.0 / (r * r * r);
      acceleration.x[i] += dx * rr3 * m[j];
      acceleration.y[i] += dy * rr3 * m[j];
      acceleration.z[i] += dz * rr3 * m[j];
      acceleration.x[j] -= dx * rr3 * m[i];
      acceleration.y[j] -= dy * rr3 * m[i];
      acceleration.z[j] -= dz * rr3 * m[i];
    };

    if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
      auto const &ch_px = particles.chain_pos().x;
      auto const &ch_py = particles.chain_pos().y;
      auto const &ch_pz = particles.chain_pos().z;
      auto const &idx = particles.index();

      size_t size = particles.number();
      for (size_t i = 0; i < size - 1; ++i)
        force(ch_px[i], ch_py[i], ch_pz[i], idx[i], idx[i + 1]);

      for (size_t i = 0; i < size - 2; ++i)
        force(ch_px[i] + ch_px[i + 1], ch_py[i] + ch_py[i + 1], ch_pz[i] + ch_pz[i + 1], idx[i], idx[i + 2]);

      for (size_t i = 0; i < size; ++i)
        for (size_t j = i + 3; j < size; ++j)
          force(px[idx[j]] - px[idx[i]], py[idx[j]] - py[idx[i]], pz[idx[j]] - pz[idx[i]], idx[i], idx[j]);

    } else {
      for (size_t i = 0; i < num; ++i) {
        for (size_t j = i + 1; j < num; ++j) {
          force(px[j] - px[i], py[j] - py[i], pz[j] - pz[i], i, j);
        }
      }
    }

  }
}
#endif //SPACEHUB_NEWTONIAN_HPP
