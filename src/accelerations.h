//
// Created by yihan on 2/26/19.
//

#ifndef SPACEHUB_ACCELERATIONS_H
#define SPACEHUB_ACCELERATIONS_H

#include "dev-tools.h"
#include "core-computation.h"

namespace SpaceH {

    template<typename Coord, bool IsVelDep>
    class Accelerations {
    public:

        SPACEHUB_STD_SCALAR_INTERFACES(acc, acc_);

        Accelerations() = default;

        explicit Accelerations(size_t size) {
            acc_.x.resize(size);
            acc_.y.resize(size);
            acc_.z.resize(size);
        };

        void resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, acc_.x, acc_.y, acc_.z);
        }

        void reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, acc_.x, acc_.y, acc_.z);
        }

        size_t number() {
            return acc_.x.size();
        }

    private:
        Coord acc_;
    };


    template<typename Coord>
    class Accelerations<Coord, true> {
    public:

        SPACEHUB_STD_SCALAR_INTERFACES(acc, acc_);
        SPACEHUB_STD_SCALAR_INTERFACES(v_dep_acc, v_dep_acc_);
        SPACEHUB_STD_SCALAR_INTERFACES(v_indep_acc, v_indep_acc_);

        Accelerations() = default;

        explicit Accelerations(size_t size)
            : ax_(size), ay_(size), az_(size), vd_ax_(size), vd_ay_(size), vd_az_(size), vid_ax_(size), vid_ay_(size), vid_az_(size) {}

        void resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, ax_, ay_, az_, vd_ax_, vd_ay_, vd_az_, vid_ax_, vid_ay_, vid_az_);
        }

        void reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, ax_, ay_, az_, vd_ax_, vd_ay_, vd_az_, vid_ax_, vid_ay_, vid_az_);
        }

        size_t number() {
            return ax_.size();
        }

    private:
        Coord acc_;
        Coord v_dep_acc_;
        Coord v_indep_acc_;
    };
}
#endif //SPACEHUB_ACCELERATIONS_H
