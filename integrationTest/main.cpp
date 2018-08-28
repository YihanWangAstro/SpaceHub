#include "../spaceHub.h"
#include <iomanip>
using namespace SpaceH;
//using scalar = double;
using scalar = precise_d;
const size_t N = 2;//SpaceH::DYNAMICAL;
int main(int argc, char **argv) {

    //using f = KarmackNewtonian<scalar,N>;
    using force = NewtonianForce<scalar, N>;
    //using PN = PostNewtonianForce<scalar,N,true,false,false>;

    //using sys = Basic<force>;
    using sys = ARchain<force>;

    using simulation = Nbody<sys>;
    //using simulation = Nbody<sys, IAS15, GaussDadau>;
    //using simulation = Nbody<sys, constIterator, GaussDadau>;
    //using simulation = Nbody<sys, constIterator, symplectic2th>;

    simulation nbody;

    //nbody.loadText("solar_earth.init");
    nbody.loadText("circular.init");
    //nbody.loadText("elliptic.init");
    //nbody.loadText("Kozai.init");

    simulation::RunArgs args;
    args.endTime = 50000 * Unit::YEAR;
    //args.registerPreOption(CallBack::DefaultWriter<sys>("out.dat", args.endTime));
    args.registerPostOption(CallBack::ShowProgressBar<sys>(args.endTime));

    nbody.run(args);
    return 0;
}
