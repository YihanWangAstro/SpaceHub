

#include "../spaceHub.h"
#include <iostream>




/*public RunArgs
{
    std::string prefix{"text_"};
    std::string path{""};
    double end_time{0};
    double out_dens{1000};
    double r_tol{0};
    double a_tol{0};
    double time_step{0};
}

template<typename dtype, size_t N, template<typename,size_t>class Force>
void SimpleTest(RunArgs& args)
{
    using Grav    = Force<dtype,N>;
    using archain = SpaceH::ARchain<Force>;
    using GAR     = SpaceH::GAR<Force>;
    using simple  = SpaceH::SimpleSys<Force>;
    
    SpaceH::Nbody<archain> nbody1;
    SpaceH::Nbody<GAR>     nbody2;
    SpaceH::Nbody<simple>  nbody3;
    
    nbody1.loadText(args.path.c_str());
    nbody2.loadText(args.path.c_str());
    nbody3.loadText(args.path.c_str());
    
    nbody1.
    
}*/
using scalar = double;//SpaceH::precise_d;
//using scalar = SpaceH::precise_d;
const size_t N = 2;//SpaceH::DYNAMICAL;
int main(int argc, char**argv)
{
    //std::cout << std::scientific << std::setprecision(16);// << nbody.particles;
    
    //using f = SpaceH::KarmackNewtonian<scalar,N>;
    using force = SpaceH::NewtonForce<scalar,N>;
    using PNf = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    
    using sys = SpaceH::Basic<force>;
    //using sys = SpaceH::ARchain<force>;
    //SpaceH::Nbody<sys> nbody;
    //SpaceH::Nbody<sys, SpaceH::IAS15, SpaceH::GaussDadau> nbody;
    SpaceH::Nbody<sys, SpaceH::constIterator, SpaceH::symplectic2th> nbody;
    
    
    
    //nbody.loadText("solar_earth.init");
    nbody.loadText("circular.init");
  
    
   // nbody.iterator.setRelativeError(1e-6);

    double timeLimit = 1000000*SpaceH::Unit::YEAR;
    double dt = 0.01*SpaceH::Unit::YEAR;
    double tout = 0;

    std::ofstream os("out.dat");
    std::ofstream eos("out.err");
    os.sync_with_stdio(false);
    eos.sync_with_stdio(false);
    os  << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);
    std::cout << std::scientific << std::setprecision(16);
   
    SpaceH::Timer clock;
    
    clock.start();
    
    for(; nbody.particles.time() < timeLimit;)
    {
        nbody.advanceOneStep();
       // std::cout << nbody.particles.time() << '\n';
        /*if(nbody.particles.time() > tout)
        {
            os  << nbody.particles;
            eos << nbody.particles.time()/SpaceH::Unit::YEAR << ' '
                << getTotalEnergy(nbody.particles.mass(), nbody.particles.pos(), nbody.particles.vel()) << ' '
                //<< nbody.particles.omega() << ' '
                //<< nbody.particles.bindE() << ' '
                //<< nbody.stepLength <<' '
                //<< SpaceH::getEnergyErr(nbody.particles.mass(),nbody.particles.pos(),nbody.particles.vel(),nbody.particles.bindE()) <<"\r\n";
                << "\r\n";
            tout += dt;
        }*/
    }
    //std::cout << nbody.steps << '\n';
    
    std::cout << std::setprecision(6) << clock.getTime() << std::endl;
    os.close();
    eos.close();
    return 0;
}
