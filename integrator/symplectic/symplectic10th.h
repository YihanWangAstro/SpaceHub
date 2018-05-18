////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_10TH_INTEGRATOR_H
#define SYMPLECTIC_10TH_INTEGRATOR_H
template <typename ParticSys>
class symplectic10th
{
    //////////////////////////////////Interface/////////////////////////////////////
public:
    static const int order{10};
    void integrate(ParticSys& particles, double stepLength);
};
    /////////////////////////////Implement Function/////////////////////////////////
template <typename ParticSys>
void symplectic10th<ParticSys>::integrate(ParticSys& particles, double stepLength)
{
    particles.AdvancePos(0.5*stepLength);
    particles.AdvanceVel(stepLength);
    particles.AdvancePos(0.5*stepLength);
}
#endif
