
#ifndef PARTICLES_H
#define PARTICLES_H

#include "macros.h"
#include "type_class.h"
#include "dev_tools.h"
#include "core_computation.h"

namespace SpaceH {

    template<typename TypeClass>
    class BasicParticle {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeClass);
        /* Typedef */
        ScalarArray px;
        ScalarArray py;
        ScalarArray pz;
        ScalarArray vx;
        ScalarArray vy;
        ScalarArray vz;
        ScalarArray mass;
        IndexArray  idn;
    };

    template<typename TypeClass>
    class FiniteSizeParticles {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeClass);
        /* Typedef */
        ScalarArray px;
        ScalarArray py;
        ScalarArray pz;
        ScalarArray vx;
        ScalarArray vy;
        ScalarArray vz;
        ScalarArray mass;
        ScalarArray radius;
        IndexArray  idn;
    };
}
#endif

