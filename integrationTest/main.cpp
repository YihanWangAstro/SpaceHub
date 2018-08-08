#include "../wrapper/singleRun.h"

#include "../wrapper/errCurve.h"
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
const size_t N = SpaceH::DYNAMICAL;
int main(int argc, char**argv)
{
    //std::cout << std::scientific << std::setprecision(16);// << nbody.particles;
    
    
    //using f = SpaceH::KarmackNewtonian<scalar,N>;
    /*using force = SpaceH::NewtonForce<scalar,N>;
    using PNf = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    
    using sys = SpaceH::ARchain<force,PNf>;
    
   SpaceH::Nbody<sys> nbody;
    
    
    
    nbody.loadText("kepler2.init");
  
    nbody.setStepLength(0.0001*YEAR);
    
   // nbody.iterator.setRelativeError(1e-6);

    double timeLimit = 10*YEAR;
    double dt = 0.001*YEAR;
    double tout = 0;

    std::ofstream os("out.dat");
    std::ofstream eos("out.eng");
    os.sync_with_stdio(false);
    eos.sync_with_stdio(false);
    os  << std::scientific << std::setprecision(16);
    eos << std::scientific << std::setprecision(16);
    
   
    SpaceH::Timer clock;
    
    clock.start();
    
    for(; nbody.particles.time() < timeLimit;)
    {
        nbody.advanceOneStep();
        //std::cout << nbody.particles.time << '\n';
        if(nbody.particles.time() > tout)
        {
            os  << nbody.particles;
            eos << nbody.particles.time()/YEAR << ' '
                << getTotalEnergy(nbody.particles.mass(), nbody.particles.pos(), nbody.particles.vel()) << ' '
                << nbody.particles.omega() << ' '
                << nbody.particles.bindE() << ' '
                << nbody.stepLength <<' '
                << SpaceH::getEnergyErr(nbody.particles.mass(),nbody.particles.pos(),nbody.particles.vel(),nbody.particles.bindE()) <<"\r\n";
            tout += dt;
        }
    }
    //std::cout << nbody.steps << '\n';
    
    std::cout << std::setprecision(6) << clock.getTime() << std::endl;
    os.close();
    eos.close();
    return 0;*/
    using f = SpaceH::NewtonForce<scalar,N>;
    using pn = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    std::cout << std::scientific << std::setprecision(6);
    
    SpaceH::RunArgs arg;
    std::string prefix = "NO";
    {
        arg.outputDensity = 20000;
        
        arg.endTime = 2000*YEAR;
        
        arg.path = "mercury.init";
        arg.name = prefix + "mercury";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = "kepler0.init";
        arg.name = prefix + "kepler0";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = "kepler.init";
        arg.name = prefix + "kepler";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.endTime = 20*YEAR;
        arg.path = "kepler2.init";
        arg.name = prefix + "kepler2";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
    }
    
    {
        arg.endTime = 2000*YEAR;
        
        arg.path = "mercury.init";
        arg.name = prefix + "Gmercury";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = "kepler0.init";
        arg.name = prefix + "Gkepler0";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = "kepler.init";
        arg.name = prefix + "Gkepler";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.endTime = 20*YEAR;
        arg.path = "kepler2.init";
        arg.name = prefix + "Gkepler2";
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
    }
    
    {
        arg.outputDensity = 20000;
        
        arg.endTime = 2000*YEAR;
        
        arg.path = "mercury.init";
        arg.name = prefix + "mercury";
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = "kepler0.init";
        arg.name = prefix + "kepler0";
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = "kepler.init";
        arg.name = prefix + "kepler";
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.endTime = 20*YEAR;
        arg.path = "kepler2.init";
        arg.name = prefix + "kepler2";
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
    }
    
    {
        arg.endTime = 2000*YEAR;
        
        arg.path = "mercury.init";
        arg.name = prefix + "Gmercury";
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = "kepler0.init";
        arg.name = prefix + "Gkepler0";
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = "kepler.init";
        arg.name = prefix + "Gkepler";
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.endTime = 20*YEAR;
        arg.path = "kepler2.init";
        arg.name = prefix + "Gkepler2";
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
    }
}
