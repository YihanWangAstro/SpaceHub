//
// Created by yihan on 2/26/19.
//

#ifndef SPACEHUB_ACCELERATIONS_H
#define SPACEHUB_ACCELERATIONS_H

#include "dev-tools.h"
#include "core-computation.h"

namespace SpaceH {

    template<typename ScalarArray, bool IsVelDep>
    class Accelerations {
    public:

        SPACEHUB_STD_ARRAY_INTERFACES(ax, ax_);
        SPACEHUB_STD_ARRAY_INTERFACES(ay, ay_);
        SPACEHUB_STD_ARRAY_INTERFACES(az, az_);

        Accelerations() = default;

        explicit Accelerations(size_t size) : ax_(size), ay_(size), az_(size) {};

        void resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, ax_, ay_, az_);
        }

        void reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, ax_, ay_, az_);
        }

        size_t number() {
            return ax_.size();
        }

    private:
        ScalarArray ax_;
        ScalarArray ay_;
        ScalarArray az_;
    };


    template<typename ScalarArray>
    class Accelerations<ScalarArray, true> {
    public:

        SPACEHUB_STD_ARRAY_INTERFACES(ax, ax_);
        SPACEHUB_STD_ARRAY_INTERFACES(ay, ay_);
        SPACEHUB_STD_ARRAY_INTERFACES(az, az_);
        SPACEHUB_STD_ARRAY_INTERFACES(vd_ax, vd_ax_);
        SPACEHUB_STD_ARRAY_INTERFACES(vd_ay, vd_ay_);
        SPACEHUB_STD_ARRAY_INTERFACES(vd_az, vd_az_);
        SPACEHUB_STD_ARRAY_INTERFACES(vid_ax, vid_ax_);
        SPACEHUB_STD_ARRAY_INTERFACES(vid_ay, vid_ay_);
        SPACEHUB_STD_ARRAY_INTERFACES(vid_az, vid_az_);

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
        ScalarArray ax_;
        ScalarArray ay_;
        ScalarArray az_;
        ScalarArray vd_ax_;
        ScalarArray vd_ay_;
        ScalarArray vd_az_;
        ScalarArray vid_ax_;
        ScalarArray vid_ay_;
        ScalarArray vid_az_;
    };


    template<typename Acc>
    void sum_all_acc(Acc &acc) {
        calc::array_add(acc.ax(), acc.vid_ax(), acc.vd_ax());
        calc::array_add(acc.ay(), acc.vid_ay(), acc.vd_ay());
        calc::array_add(acc.az(), acc.vid_az(), acc.vd_az());
    }
}
#endif //SPACEHUB_ACCELERATIONS_H
