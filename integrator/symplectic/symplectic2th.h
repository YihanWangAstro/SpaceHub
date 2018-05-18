////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:size_tegrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_2TH_INTEGRATOR_H
#define SYMPLECTIC_2TH_INTEGRATOR_H
template <typename ParticSys>
class symplectic2th
{
    //////////////////////////////////size_terface/////////////////////////////////////
    typedef typename ParticSys::Scalar Scalar;
public:
    static const int order{2};
    void integrate(ParticSys& particles, Scalar stepLength);
};
    /////////////////////////////Implement Function/////////////////////////////////
template <typename ParticSys>
void symplectic2th<ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
{
    particles.advancePos(0.5*stepLength);
    particles.advanceVel(stepLength);
    particles.advancePos(0.5*stepLength);
}
#endif
