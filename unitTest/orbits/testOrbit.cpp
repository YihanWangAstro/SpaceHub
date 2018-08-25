#include <fstream>
#include "../../initCreator/orbits.h"
#include "../../initCreator/InitFileCreator.h"

using Scalar    = double;
using particle  = SpaceH::Particle<Scalar>;
using Creator   = SpaceH::InitCreator<particle>;
using KObt = SpaceH::Kepler<Scalar>;

constexpr size_t N = 5000;
constexpr size_t orbitsNum = 10;
constexpr size_t interval = N / orbitsNum;

int main() {
    {
        Scalar e = 0.6;
        Scalar p = fabs(1 - e * e) * SpaceH::Unit::AU;;

        KObt orbit(1, 1e-3, p, e);
        Creator crt;

        for (size_t i = 0; i < N; ++i) {
            if (i % interval == 0)
                orbit.randomAngles();

            //orbit.randomPhase(-SpaceH::Unit::YEAR, 0);
            orbit.randomPhase();

            crt.addParticle(orbit.primary());
            crt.addParticle(orbit.secondary());
        }
        crt.writeToFile("elliptic.init");
    }

    {
        Scalar e = 1;
        Scalar p = 0.5 * SpaceH::Unit::AU;;

        KObt orbit(1, 1e-3, p, e);
        Creator crt;

        for (size_t i = 0; i < N; ++i) {
            if (i % interval == 0)
                orbit.randomAngles();

            orbit.randomPhase(-SpaceH::Unit::YEAR, SpaceH::Unit::YEAR);
            //orbit.randomPhase();

            crt.addParticle(orbit.primary());
            crt.addParticle(orbit.secondary());
        }
        crt.writeToFile("parabolic.init");
    }

    {
        Scalar e = 2.5;
        Scalar p = fabs(1 - e * e) * SpaceH::Unit::AU;;

        KObt orbit(1, 1e-3, p, e);
        Creator crt;

        for (size_t i = 0; i < N; ++i) {
            if (i % interval == 0)
                orbit.randomAngles();

            orbit.randomPhase(-SpaceH::Unit::YEAR, SpaceH::Unit::YEAR);
            // orbit.randomPhase();

            crt.addParticle(orbit.primary());
            crt.addParticle(orbit.secondary());
        }
        crt.writeToFile("hyperbolic.init");
    }
}
