////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_8TH_INTEGRATOR_H
#define SYMPLECTIC_8TH_INTEGRATOR_H
template <typename ParticSys>
class symplectic8th
{
    //////////////////////////////////Interface/////////////////////////////////////
public:
    static const int order{8};
    void integrate(ParticSys& particles, double stepLength);
};
    /////////////////////////////Implement Function/////////////////////////////////
template <typename ParticSys>
void symplectic8th<ParticSys>::integrate(ParticSys& particles, double stepLength)
{
    particles.AdvancePos(0.5*stepLength);
    particles.AdvanceVel(stepLength);
    particles.AdvancePos(0.5*stepLength);
}
#endif
