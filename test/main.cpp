#include "../spaceX.h"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std::chrono;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> resolutionClock;

//using scalar = float;//SpaceH::kahan<double>;
using scalar = SpaceH::kahan<double>;

const size_t N = 3;
int main(int argc, char**argv)
{
    //std::cout << std::scientific << std::setprecision(16);// << nbody.particles;
    
    resolutionClock start;
    resolutionClock now;
    
    
    using type = SpaceH::ProtoType<scalar,N>;
    
    //using f = SpaceH::KarmackNewtonian<scalar,N>;
    using f = SpaceH::NewtonForce<scalar,N>;
    using PNf = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    using force = SpaceH::VelIndepForce< f,scalar,N>;
    using PNforce = SpaceH::VelDepForce< PNf,scalar,N>;
    using no = SpaceH::EmptyForce<scalar,N>;
    using interaction = SpaceH::Interaction<force,no,no,no>;
    
    using particle = ChainParticles<scalar,N,interaction::isVelDep>;
    using regu = logH<particle>;
    using ps = ReguSystem<particle, interaction,regu>;
    
    //using particle = Particles<scalar,N,interaction::isVelDep>;
    //using ps = particleSystem<particle, interaction>;
    
    using integ = symplectic2th<ps>;
    using it = BSIterator<ps,integ>;
    //using it = constIterator<ps,integ>;
    
    using ds = dynamicSystem<ps,it>;
    
    
    ds nbody;
    
    nbody.loadText("kepler2.init");
  
    nbody.setStepLength(0.0001*YEAR);
    
   // nbody.iterator.setRelativeError(1e-6);

    double timeLimit = 50*YEAR;
    double dt = 5*YEAR;
    double tout = 0;

    /*std::ofstream os("out.dat");
    std::ofstream eos("out.eng");
    os.sync_with_stdio(false);
    eos.sync_with_stdio(false);
    os  << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);*/
    
    /*typename type::ScalarBuffer mapp;
    mapp << nbody.particles;
    SpaceH::print(mapp);
    std::cout << mapp.size() << '\n';*/
    //std::cout <<nbody.particles;
    
   // nbody.advanceOneStep();
    
    //std::cout <<nbody.particles;
    //std::cout << nbody.particles << std::endl;
    /*std::cout << nbody.particles.dynState.volume() << std::endl;
    for(int i = 0 ; i < nbody.particles.dynState.volume(); i++)
    {
        std::cout << nbody.particles.array()[i] << ' ';
    }*/
    
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
                 << 1 << ' ' << 1 << ' ' << 1 << ' '
                //<< nbody.particles.omega() << ' '
                //<< nbody.particles.bindE() << ' '
            << nbody.stepLength <<' ' << "\r\n";
                //<< SpaceH::getEnergyErr(nbody.particles.mass(),nbody.particles.pos(),nbody.particles.vel(),nbody.particles.bindE()) <<"\r\n";
                //<< nbody.iterator.iterDepth << "\r\n";
            tout += dt;
        }*/
    }
    //std::cout << nbody.steps << '\n';
    now           = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(now - start);
    auto runtime  = double(duration.count()) * milliseconds::period::num / milliseconds::period::den;
    std::cout << std::setprecision(6) << runtime << std::endl;
    //os.close();
    //eos.close();
    return 0;
}
