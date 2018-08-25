#include "../spaceHub.h"
#include <iomanip>
//using scalar = double;
using scalar = SpaceH::precise_d;
const size_t N = 3;//SpaceH::DYNAMICAL;
int main(int argc, char **argv) {

    //using f = SpaceH::KarmackNewtonian<scalar,N>;
    using force = SpaceH::NewtonForce<scalar, N>;
    //using PN = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;

    //using sys = SpaceH::Basic<force>;
    using sys = SpaceH::ARchain<force>;
    using simulation = SpaceH::Nbody<sys>;

    //using simulation = SpaceH::Nbody<sys, SpaceH::IAS15, SpaceH::GaussDadau>;
    //using simulation = SpaceH::Nbody<sys, SpaceH::constIterator, SpaceH::GaussDadau>;
    //using simulation = SpaceH::Nbody<sys, SpaceH::constIterator, SpaceH::symplectic2th>;

    simulation nbody;

    //nbody.setStepLength(0.01*SpaceH::Unit::YEAR);
    //nbody.loadText("solar_earth.init");
    //nbody.loadText("circular.init");
    //nbody.loadText("elliptic.init");
    nbody.loadText("Kozai.init");

    //nbody.iterator.setRelativeError(1e-6);

    simulation::RunArgs args;
    args.endTime = 500 * SpaceH::Unit::YEAR;


    std::ofstream os("out.dat");
    std::ofstream eos("out.err");
    os << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);
    std::ios_base::sync_with_stdio(false);


    scalar dt = 0.01 * SpaceH::Unit::YEAR;
    scalar tout = 0;

    auto output = [&](sys &partc) {
        if (partc.time() > tout) {
            os << partc;
            eos << partc.time() / SpaceH::Unit::YEAR << ' ' << partc.totalEnergy() << "\r\n";
            tout += dt;
        }
    };

    args.registerPreOption(output);

    SpaceH::Progress bar(args.endTime);

    args.registerPostOption([&](sys &partc) { bar.autoShow(partc.time()); });

    nbody.run(args);

    os.close();
    eos.close();
    return 0;
}
