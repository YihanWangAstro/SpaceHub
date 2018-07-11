#include "../spaceX.h"
#include <chrono>
#include <fstream>
#include <iostream>
//bool NOTICE::Message = true;
using namespace std::chrono;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> resolutionClock;

using scalar = double;
int main(int argc, char**argv)
{
    //std::cout << std::scientific << std::setprecision(16);// << nbody.particles;
    resolutionClock start;
    resolutionClock now;
    start         = high_resolution_clock::now();
    //spaceX<VelDepSystem<2,PN1th,TTL>, symplectic2th, BSIterator> nbody;
    //spaceX<NewtonianSystem<2>, symplectic2th, BSIterator> nbody;
    //using type  = typeSet::type<double,2>;
    
    using particle = ReguParticles<scalar,2>;
    using force = interact::NewtonGrav<scalar,2>;
    using no = interact::Empty<scalar,2>;
    using interaction = interact::Interaction<force,no,no,no>;
    
    using regu = logH<particle>;
    using ps = ReguSystem<particle, interaction, regu>;
    
    using integ = symplectic2th<ps>;
    using it = BSIterator<ps,integ>;
    //using it = constIterator<ps,integ>;
    
    using ds = dynamicSystem<ps,it>;
    
    
    ds nbody;
    
    nbody.loadText("kepler.init");
  
    nbody.setStepLength(0.01*YEAR);
    
   // nbody.iterator.setRelativeError(1e-6);

    double timeLimit = 1000*YEAR;
    double dt = 1*YEAR;
    double tout = 0;

    std::ofstream os("out.dat");
    std::ofstream eos("out.eng");
    os.sync_with_stdio(false);
    eos.sync_with_stdio(false);
    os  << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);
    
    
    //std::cout << nbody.particles << std::endl;
    /*std::cout << nbody.particles.dynState.volume() << std::endl;
    for(int i = 0 ; i < nbody.particles.dynState.volume(); i++)
    {
        std::cout << nbody.particles.array()[i] << ' ';
    }*/
    
    //nbody.advanceOneStep();
    for(; nbody.particles.time() < timeLimit;)
    {
        nbody.advanceOneStep();
        //std::cout << nbody.particles.time << '\n';
        /*if(nbody.particles.time() > tout)
        {
            os  << nbody.particles;
            eos << nbody.particles.time()/YEAR << ' '
                << getTotalEnergy(nbody.particles.mass(), nbody.particles.pos(), nbody.particles.vel()) << ' '
                << getPotentialEnergy(nbody.particles.mass(), nbody.particles.pos()) << ' '
                << getKineticEnergy(nbody.particles.mass(), nbody.particles.vel()) << ' '
            << nbody.stepLength << "\r\n";
                //<< nbody.iterator.iterDepth << "\r\n";
            tout += dt;
        }*/
    }
   // std::cout << nbody.steps << '\n';
    now           = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(now - start);
    auto runtime  = double(duration.count()) * milliseconds::period::num / milliseconds::period::den;
    std::cout << std::setprecision(6) << runtime << std::endl;
    os.close();
    eos.close();
    return 0;
}
