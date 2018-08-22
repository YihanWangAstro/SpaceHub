#include "../spaceHub.h"
#include <iostream>

//using scalar = double;
using scalar = SpaceH::precise_d;
const size_t N = 3;//SpaceH::DYNAMICAL;
int main(int argc, char**argv)
{
    
    //using f = SpaceH::KarmackNewtonian<scalar,N>;
    using force = SpaceH::NewtonForce<scalar,N>;
    //using PN = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    
    //using sys = SpaceH::Basic<force>;
    using sys = SpaceH::ARchain<force>;
    SpaceH::Nbody<sys> nbody;
    //SpaceH::Nbody<sys, SpaceH::IAS15, SpaceH::GaussDadau> nbody;
    //SpaceH::Nbody<sys, SpaceH::constIterator, SpaceH::GaussDadau> nbody;
    //SpaceH::Nbody<sys, SpaceH::constIterator, SpaceH::symplectic2th> nbody;
    
    //nbody.setStepLength(0.01*SpaceH::Unit::YEAR);
    //nbody.loadText("solar_earth.init");
    //nbody.loadText("circular.init");
    //nbody.loadText("elliptic.init");
      nbody.loadText("Kozai.init");
    
    //nbody.iterator.setRelativeError(1e-6);

    scalar timeLimit = 100*SpaceH::Unit::YEAR;
    scalar dt = 0.0005*SpaceH::Unit::YEAR;
    scalar tout = 0;

    std::ofstream os("out.dat");
    std::ofstream eos("out.err");

    os  << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);
    std::cout << std::scientific << std::setprecision(16);
    std::ios_base::sync_with_stdio(false);
    
    SpaceH::Timer clock;
    
    clock.start();
    for(; nbody.particles.time() < timeLimit ;)
    {
        nbody.advanceOneStep();
        if(nbody.particles.time() > tout)
        {
            os  << nbody.particles;
            eos << nbody.particles.time()/SpaceH::Unit::YEAR << ' '
                << nbody.particles.totalEnergy() << "\r\n";
            tout += dt;
        }
    }
    
    std::cout << std::setprecision(6) << clock.getTime() << std::endl;
    os.close();
    eos.close();
    return 0;
}
