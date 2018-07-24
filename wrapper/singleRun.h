
#ifndef WRAPPER_SINGLE_H
#define WRAPPER_SINGLE_H


#include "../spaceHub.h"
namespace SpaceH
{
    template
    <
        typename dtype,
        size_t N,
        template<typename,size_t>class Force,
        template <typename, typename, typename, typename, template<typename> class> Method
    >
    void SimpleTest(RunArgs& args)
    {
        using Grav    = Force<dtype,N>;
        using System  = Method<Force>;
        
        SpaceH::Nbody<System> nbody;
        
        nbody1.loadText(args.path.c_str());
        
        dtype end_time = args.endTime;
        dtype output_time = end_time / args.outputDensity;
        dtype output_point = 0;
        
        std::ofstream out_file(args.prefix + args.runName + ".dat");
        std::ofstream err_file(args.prefix + args.runName + ".err");
        std::ofstream log_file(args.prefix + args.runName + ".log");
        
        SpaceH::Timer clock;
        
        clock.start();
        
        for( ;nbody.particles.time() < end_time; )
        {
            nbody.advanceOneStep();
            if(nbody.particles.time() > output_point)
            {
                out_file << nbody.particles;
                err_file << nbody.particles.time()/YEAR << ' '
                         << SpaceH::getEnergyErr(nbody.particles.mass(),nbody.particles.pos(),nbody.particles.vel(),nbody.particles.bindE())
                         << "\r\n";
                output_point += out_time;
            }
        }
        
        log_file << clock.getTime();
        
        out_file.close();
        err_file.close();
        log_file.close();
    }
}
#endif
