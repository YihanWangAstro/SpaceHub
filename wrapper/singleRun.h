
#ifndef WRAPPER_SINGLE_H
#define WRAPPER_SINGLE_H


#include "../spaceHub.h"
#include <fstream>
#include "proto.h"
namespace SpaceH
{
    template
    <
        typename BasicF,
        typename VelForce = void,
        typename ExtPosForce = void,
        typename ExtVelForce = void,
        template
        <
            typename,
            typename = void,
            typename = void,
            typename = void,
            template<typename> class Regularitor = SpaceH::LogH
        > class Algorithm = ARchain,
        template<typename, typename> class ODEiterator = SpaceH::BSIterator,
        template<typename> class           Integrator  = SpaceH::symplectic2th
    >
    void SimpleRun(RunArgs& args)
    {

        using System  = Algorithm<BasicF, VelForce, ExtPosForce, ExtVelForce>;
        using dtype   = typename BasicF::type::Scalar;
        
        SpaceH::Nbody<System, ODEiterator, Integrator> nbody;
        
        nbody.loadText(args.path.c_str());
        
        //nbody.setStepLength(0.0001*YEAR);
        dtype end_time = args.endTime;
        dtype output_time = end_time / args.outputDensity;
        dtype output_point = 0;
        
        std::ofstream out_file(args.prefix + args.name + ".dat");
        std::ofstream err_file(args.prefix + args.name + ".err");
        std::ofstream log_file(args.prefix + args.name + ".log");
        
        out_file << std::scientific << std::setprecision(16);
        err_file << std::scientific << std::setprecision(16);
        log_file << std::scientific << std::setprecision(16);
        
        SpaceH::Progress bar(end_time);
        
        bar.start();
        dtype E0 = SpaceH::getTotalEnergy(nbody.particles.mass(), nbody.particles.pos(), nbody.particles.vel());
        for( ;nbody.particles.time() < end_time; )
        {
            nbody.advanceOneStep();
            
            if(nbody.particles.time() > output_point)
            {
                out_file << nbody.particles;
                err_file << nbody.particles.time()/SpaceH::Unit::YEAR << ' '
                         //<< SpaceH::getEnergyErr(nbody.particles.mass(),nbody.particles.pos(),nbody.particles.vel(),nbody.particles.bindE())
                << SpaceH::abs((SpaceH::getTotalEnergy(nbody.particles.mass(), nbody.particles.pos(), nbody.particles.vel()) - E0)/E0)
                << ' ' << nbody.iterator.getRejRate() 
                         << "\r\n";
                output_point += output_time;
            }
            bar.autoShow(static_cast<float>(nbody.particles.time()));
        }
        
        log_file << bar.getTime();
        
        out_file.close();
        err_file.close();
        log_file.close();
    }
}
#endif
