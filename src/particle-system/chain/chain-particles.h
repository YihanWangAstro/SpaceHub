//
// Created by yihan on 3/1/19.
//

#ifndef SPACEHUB_CHAIN_PARTICLES_H
#define SPACEHUB_CHAIN_PARTICLES_H

#include "particle.h"
namespace SpaceH {

    template <typename Derived, bool IsVelDep>
    class SoAChainParticle : SoAParticles<Derived, IsVelDep>{
    public:
        DECLARE_CRTP_ACCESSOR(chain_pos, typename Derived::Coord, Derived);
        DECLARE_CRTP_ACCESSOR(chain_vel, typename Derived::Coord, Derived);

    private:
        SoAChainParticle() = default;
        friend Derived;
    };


    template <typename Derived>
    class SoAChainParticle<Derived, true> : public SoAChainParticle<Derived, false>{
    public:
        DECLARE_CRTP_ACCESSOR(chain_aux_vel, typename Derived::Coord, Derived);
    private:
        SoAChainParticle() = default;
        friend Derived;
    };
}
#endif //SPACEHUB_CHAIN_PARTICLES_H
