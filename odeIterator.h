////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:odeIterator.h                                                                                              //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ITERATOR_H
#define ITERATOR_H
template <typename Derived, typename ParticSys, typename Integrator>
class odeIterator
{
public:
    //////////////////////////////////size_terface/////////////////////////////////////
    typedef typename ParticSys::Scalar Scalar;
    Scalar iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength)
    {
        return static_cast<Derived*>(this)->impl_iterate(particles, integrator, stepLength);
    }
    void setRelativeError(Scalar relError){ relativeError = relError; }
    void setAbsoluteError(Scalar absError){ absoluteError = absError; }
    virtual ~odeIterator(){}
    ///////////////////////////////Member variables/////////////////////////////////
protected:
    Scalar absoluteError{1e-15};
    Scalar relativeError{1e-15};
};

#endif
