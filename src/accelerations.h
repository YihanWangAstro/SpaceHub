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

        SPACEHUB_STD_INTERFACES(acc, acc_);

        Accelerations() = default;

        explicit Accelerations(size_t size) : acc_(size) {};

        void resize(size_t new_sz) {
            acc_.resize(new_sz);
        }

        void reserve(size_t new_cap) {
            acc_.reserve(new_cap);
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

        SPACEHUB_STD_INTERFACES(acc, acc_);
        SPACEHUB_STD_INTERFACES(vel_dep_acc, vel_dep_acc_);
        SPACEHUB_STD_INTERFACES(vel_indep_acc, vel_indep_acc_);

        Accelerations() = default;

        explicit Accelerations(size_t size)
            : acc_(size), vel_dep_acc_(size), vel_indep_acc_(size) {}

        void resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, acc_, vel_dep_acc_, vel_indep_acc_);
        }

        void reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, acc_, vel_dep_acc_, vel_indep_acc_);
        }

        size_t number() {
            return acc_.x.size();
        }

    private:
        Coord acc_;
        Coord vel_dep_acc_;
        Coord vel_indep_acc_;
    };
}
#endif //SPACEHUB_ACCELERATIONS_H
