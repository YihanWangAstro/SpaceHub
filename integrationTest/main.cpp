#include "../spaceHub.h"
#include <iomanip>
using namespace SpaceH;
//using scalar = double;
using scalar = precise_d;
const size_t N = 2;//SpaceH::DYNAMICAL;
using type = SpaceH::ProtoType<scalar, N>;
int main(int argc, char **argv) {

    //using f = KarmackNewtonian<type>;
    using force = NewtonianForce<type>;
    //using PN = PostNewtonianForce<type,true,false,false>;

    //using sys = Basic<force>;
    //using sys =  GAR<force>;
    using sys = ARchain <force>;
    using simulation = Nbody<sys>;
    //using simulation = Nbody<sys, IAS15, GaussDadau>;
    //using simulation = Nbody<sys, constIterator, GaussDadau>;
    //using simulation = Nbody<sys, constIterator, symplectic2th>;

    simulation nbody;

    //nbody.loadText("solar_earth.init");
    nbody.loadText("circular.init");
    //nbody.loadText("elliptic.init");
    //nbody.loadText("Kozai.init");
    std::cout << std::setprecision(16) << nbody.particles;

    simulation::RunArgs args;
    args.endTime = 1000 * Unit::YEAR;

    args.registerPreOption(CallBack::DefaultWriter<sys>("circular.dat", args.endTime));
    args.registerPostOption(CallBack::ShowProgressBar<sys>(args.endTime));

    nbody.run(args);
    return 0;
}
