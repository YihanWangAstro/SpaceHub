////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:size_tegrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef size_tEGRATOR_H
#define size_tEGRATOR_H
template <typename Derived, typename ParticSys>
class integrator
{
public:
    ///////////////////////////////////Interface////////////////////////////////////
    typedef typename ParticSys::Scalar Scalar;
     //////////////////////////////////Interface/////////////////////////////////////
    void integrate(ParticSys& particles, Scalar stepLength);
    virtual ~integrator(){}
    ///////////////////////////////Member variables/////////////////////////////////
    size_t order;
};
    /////////////////////////////Implement Function/////////////////////////////////
template <typename Derived, typename ParticSys>
inline void integrator<Derived, ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
{
    static_cast<Derived*>(this)->impl_integrate(particles, stepLength);
}

#endif
