////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:constIterator.h                                                                                            //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef CONSTITERATOR_H
#define CONSTITERATOR_H
#include "../odeIterator.h"
template <typename ParticSys, typename Integrator>
class constIterator : public odeIterator<constIterator<ParticSys, Integrator>, ParticSys, Integrator>
{
public:
    typedef typename ParticSys::Scalar Scalar;
    Scalar impl_iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength)
    {
        integrator.size_tegrate(particles, stepLength);
        return stepLength;
    }
};

#endif
