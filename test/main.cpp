#include "../spaceX.h"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std::chrono;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> resolutionClock;

using scalar = double;
const size_t N = 2;
int main(int argc, char**argv)
{
    //std::cout << std::scientific << std::setprecision(16);// << nbody.particles;
    
    resolutionClock start;
    resolutionClock now;
    
    
    //using f = SpaceH::KarmackNewtonian<scalar,N>;
    using force = SpaceH::NewtonForce<scalar,N>;
   // using PNf = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    
    using sys = SpaceH::ARchain<force>;
    
    SpaceH::Nbody<sys> nbody;
    
    nbody.loadText("kepler.init");
  
    nbody.setStepLength(0.0001*YEAR);
    
   // nbody.iterator.setRelativeError(1e-6);

    double timeLimit = 100000*YEAR;
    double dt = 0.01*YEAR;
    double tout = 0;

    std::ofstream os("out.dat");
    std::ofstream eos("out.eng");
    os.sync_with_stdio(false);
    eos.sync_with_stdio(false);
    os  << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);
    
    
    start         = high_resolution_clock::now();
    

    for(; nbody.particles.time() < timeLimit;)
    {
        nbody.advanceOneStep();
        //std::cout << nbody.particles.time << '\n';
        /*if(nbody.particles.time() > tout)
        {
            os  << nbody.particles;
            eos << nbody.particles.time()/YEAR << ' '
                << getTotalEnergy(nbody.particles.mass(), nbody.particles.pos(), nbody.particles.vel()) << ' '
                // << 1 << ' ' << 1 << ' ' << 1 << ' '
                << nbody.particles.omega() << ' '
                << nbody.particles.bindE() << ' '
            << nbody.stepLength <<' '
                << SpaceH::getEnergyErr(nbody.particles.mass(),nbody.particles.pos(),nbody.particles.vel(),nbody.particles.bindE()) <<"\r\n";
                //<< nbody.iterator.iterDepth << "\r\n";
            tout += dt;
        }*/
    }
    //std::cout << nbody.steps << '\n';
    now           = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(now - start);
    auto runtime  = double(duration.count()) * milliseconds::period::num / milliseconds::period::den;
    std::cout << std::setprecision(6) << runtime << std::endl;
    os.close();
    eos.close();
    return 0;
}
