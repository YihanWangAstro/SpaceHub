//
// Created by 王艺涵 on 4/12/19.
//

#ifndef SPACEHUB_POST_NEWTONIAN_H
#define SPACEHUB_POST_NEWTONIAN_H
#include "../interaction.h"
#include "../dev-tools.h"

namespace space::Interact {

    template <bool _1st_order, bool _2nd_order, bool _2_5th_order>
    class PostNewtonianGrav : public Interactions<PostNewtonianGrav<_1st_order, _2nd_order, _2_5th_order>> {
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

        template<typename Particles>
        void impl_eval_extra_vel_dep_acc(Particles const &partc, typename Particles::Coord &acc){
            size_t num = partc.number();
            auto &px = partc.pos().x;
            auto &py = partc.pos().y;
            auto &pz = partc.pos().z;
            auto &m = partc.mass();

            calc::set_arrays_zero(acc.x, acc.y, acc.z);

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
                for (size_t i = 0; i < size - 1; ++i){
                    if constexpr (_1st_order) {
                        force(ch_px[i], ch_py[i], ch_pz[i], idx[i], idx[i + 1]);
                    }
                    if constexpr (_2nd_order) {
                        force(ch_px[i], ch_py[i], ch_pz[i], idx[i], idx[i + 1]);
                    }
                    if constexpr (_2_5th_order) {
                        force(ch_px[i], ch_py[i], ch_pz[i], idx[i], idx[i + 1]);
                    }
                }

                for (size_t i = 0; i < size - 2; ++i){
                    if constexpr (_1st_order) {
                        force(ch_px[i] + ch_px[i + 1], ch_py[i] + ch_py[i + 1], ch_pz[i] + ch_pz[i + 1], idx[i], idx[i + 2]);
                    }
                    if constexpr (_2nd_order) {
                        force(ch_px[i] + ch_px[i + 1], ch_py[i] + ch_py[i + 1], ch_pz[i] + ch_pz[i + 1], idx[i], idx[i + 2]);
                    }
                    if constexpr (_2_5th_order) {
                        force(ch_px[i] + ch_px[i + 1], ch_py[i] + ch_py[i + 1], ch_pz[i] + ch_pz[i + 1], idx[i], idx[i + 2]);
                    }
                }

                for (size_t i = 0; i < size; ++i) {
                    for (size_t j = i + 3; j < size; ++j) {
                        if constexpr (_1st_order) {
                            force(px[idx[j]] - px[idx[i]], py[idx[j]] - py[idx[i]], pz[idx[j]] - pz[idx[i]], idx[i], idx[j]);
                        }
                        if constexpr (_2nd_order) {
                            force(px[idx[j]] - px[idx[i]], py[idx[j]] - py[idx[i]], pz[idx[j]] - pz[idx[i]], idx[i], idx[j]);
                        }
                        if constexpr (_2_5th_order) {
                            force(px[idx[j]] - px[idx[i]], py[idx[j]] - py[idx[i]], pz[idx[j]] - pz[idx[i]], idx[i], idx[j]);
                        }
                    }
                }
            } else {
                for (size_t i = 0; i < num; ++i) {
                    for (size_t j = i + 1; j < num; ++j) {
                        if constexpr (_1st_order) {
                            force(px[j] - px[i], py[j] - py[i], pz[j] - pz[i], i, j);
                        }
                        if constexpr (_2nd_order) {
                            force(px[j] - px[i], py[j] - py[i], pz[j] - pz[i], i, j);
                        }
                        if constexpr (_2_5th_order) {
                            force(px[j] - px[i], py[j] - py[i], pz[j] - pz[i], i, j);
                        }
                    }
                }
            }
        }

    private:
        template<typename Scalar>
        void add_1st_order_acc(Scalar dx, Scalar dy, Scalar dz) {

        }

        template<typename Particles>
        void add__2nd_order_acc(Particles const &partc, typename Particles::Coord &acc) {

        }

        template<typename Particles>
        void add_2_5th_order_acc(Particles const &partc, typename Particles::Coord &acc) {

        }
    };
}
#endif //SPACEHUB_POST_NEWTONIAN_H
