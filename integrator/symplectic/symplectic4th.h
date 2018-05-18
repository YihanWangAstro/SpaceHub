////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_4TH_INTEGRATOR_H
#define SYMPLECTIC_4TH_INTEGRATOR_H
template <typename ParticSys>
class symplectic4th
{
    //////////////////////////////////Interface/////////////////////////////////////
public:
   static const int order{4};
    void integrate(ParticSys& particles, double stepLength);
};
    /////////////////////////////Implement Function/////////////////////////////////
template <typename ParticSys>
void symplectic4th<ParticSys>::integrate(ParticSys& particles, double stepLength)
{
    particles.AdvancePos(6.7560359597983000E-1*stepLength);
    particles.AdvanceVel(1.3512071919596600E0*stepLength);
    particles.AdvancePos(-1.7560359597983000E-1*stepLength);
    particles.AdvanceVel(-1.7024143839193200E0*stepLength);
    particles.AdvancePos(-1.7560359597983000E-1*stepLength);
    particles.AdvanceVel(1.3512071919596600E0*stepLength);
    particles.AdvancePos(6.7560359597983000E-1*stepLength);
}
#endif
