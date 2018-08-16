#include "../wrapper/singleRun.h"

#include "../wrapper/errCurve.h"
#include <iostream>
using namespace SpaceH;
std::string files[]={"circular","elliptic","solar_earth","circular"};
using scalar = SpaceH::precise_d;
const size_t N = 2;//SpaceH::DYNAMICAL;
int main(int argc, char**argv)
{
    //std::cout << std::scientific << std::setprecision(16);// << nbody.particles;
    using f = SpaceH::NewtonForce<scalar,N>;
    using pn = SpaceH::PostNewtonianForce<scalar,N,true,false,false>;
    std::cout << std::scientific << std::setprecision(6);
    
    SpaceH::RunArgs arg;
    std::string prefix = "NO";
    {
        arg.outputDensity = 20000;
        
        arg.endTime = 2000*Unit::YEAR;
        
        arg.path = files[0] + ".init";
        arg.name = prefix + files[0];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = files[1] + ".init";
        arg.name = prefix + files[1];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = files[2] + ".init";
        arg.name = prefix + files[2];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.endTime = 20*Unit::YEAR;
        arg.path = files[3] + ".init";
        arg.name = prefix + files[3];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::ARchain>(arg);
    }
    
    {
        arg.endTime = 2000*Unit::YEAR;
        
        arg.path = files[0] + ".init";
        arg.name = prefix + "G" + files[0];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = files[1] + ".init";
        arg.name = prefix + "G" + files[1];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = files[2] + ".init";
        arg.name = prefix + "G" + files[2];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.endTime = 20*Unit::YEAR;
        arg.path = files[3] + ".init";
        arg.name = prefix + "G" + files[3];
        SpaceH::SimpleRun<f,void,void,void,SpaceH::GAR>(arg);
    }
    
    {
        arg.outputDensity = 20000;
        
        arg.endTime = 2000*Unit::YEAR;
        
        arg.path = files[0] + ".init";
        arg.name = prefix + files[0];
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = files[1] + ".init";
        arg.name = prefix + files[1];
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.path = files[2] + ".init";
        arg.name = prefix + files[2];
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
        
        arg.endTime = 20*Unit::YEAR;
        arg.path = files[3] + ".init";
        arg.name = prefix + files[3];
        SpaceH::performance<f,void,void,void,SpaceH::ARchain>(arg);
    }
    
    {
        arg.endTime = 2000*Unit::YEAR;
        
        arg.path = files[0] + ".init";
        arg.name = prefix + "G" + files[0];
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = files[1] + ".init";
        arg.name = prefix + "G" + files[1];
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.path = files[2] + ".init";
        arg.name = prefix + "G" + files[2];
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
        
        arg.endTime = 20*Unit::YEAR;
        arg.path = files[3] + ".init";
        arg.name = prefix + "G" + files[3];
        SpaceH::performance<f,void,void,void,SpaceH::GAR>(arg);
    }
}
