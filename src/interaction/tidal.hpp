//
// Created by 王艺涵 on 4/6/19.
//

#ifndef SPACEHUB_TIDAL_H
#define SPACEHUB_TIDAL_H

#include "../interaction.hpp"
#include "../dev-tools.hpp"
#include "../macros.hpp"
#include "../orbits/orbits.hpp"

namespace space {
  class TidalForce : public Interactions<TidalForce> {
  private:
    CREATE_METHOD_CHECK(chain_pos);

    CREATE_METHOD_CHECK(chain_vel);

    CREATE_METHOD_CHECK(index);
  public:

    template<typename Particles>
    void impl_eval_acc(Particles const &partc, typename Particles::Coord &acc) {
      impl_eval_newtonian_acc(partc, acc);
    }

    template<typename Particles>
    void impl_eval_extra_vel_dep_acc(Particles const &partc, typename Particles::Coord &acc) {

    }

    template<typename Particles>
    void impl_eval_newtonian_acc(Particles const &partc, typename Particles::Coord &acc) {
      size_t num = partc.number();
      auto &px = partc.pos().x;
      auto &py = partc.pos().y;
      auto &pz = partc.pos().z;
      auto &m = partc.mass();

      calc::set_arrays_zero(acc.x, acc.y, acc.z);

      auto force = [&](auto dx, auto dy, auto dz, auto i, auto j) {
        auto r = sqrt(dx * dx + dy * dy + dz * dz);
        auto rr3 = 1.0 / (r * r * r);
        acc.x[i] += dx * rr3 * m[j];
        acc.y[i] += dy * rr3 * m[j];
        acc.z[i] += dz * rr3 * m[j];
        acc.x[j] -= dx * rr3 * m[i];
        acc.y[j] -= dy * rr3 * m[i];
        acc.z[j] -= dz * rr3 * m[i];
      };

      if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, chain_vel) &&
                    HAS_METHOD(Particles, index)) {
        auto const &ch_px = partc.chain_pos().x;
        auto const &ch_py = partc.chain_pos().y;
        auto const &ch_pz = partc.chain_pos().z;
        auto const &ch_vx = partc.chain_vel().x;
        auto const &ch_vy = partc.chain_vel().y;
        auto const &ch_vz = partc.chain_vel().z;
        auto const &idx = partc.index();

        size_t size = partc.number();
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

  private:
    template<typename Scalar>
    auto tidal_T(Scalar eta) {
      auto x = log10(eta);
      constexpr Scalar A{-0.517};
      constexpr Scalar B{-0.906};
      constexpr Scalar C{23.88};
      constexpr Scalar D{-93.49};
      constexpr Scalar E{112.3};
      constexpr Scalar F{-44.15};
      auto logT = (((((F * x) + E) * x + D) * x + C) * x + B) * x + A;
      return pow(10, logT);
    }

    template<typename Scalar>
    auto calc_e_loss(Scalar mt, Scalar mp, Scalar Rt, Scalar Rp, Scalar a, Scalar e) {
      auto rp = a * (1 - e);
      auto e1 = (space::consts::G * mt * mt / Rt);
      auto e2 = mp * mp / mt;
      auto e3 = Rt / rp;
      auto eta = sqrt(mt / (mt + mp)) * pow(rp / Rt, 1.5);
      return e1 * e2 * e2 * e3 * e3 * e3 * e3 * e3 * e3 * tidal_T(eta);
    }
  };
}
#endif //SPACEHUB_TIDAL_H
