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
template <typename ParticSys, typename Integrator>
class constIterator
{
public:
    typedef typename ParticSys::Scalar Scalar;
    Scalar iterate(ParticSys& particles, Integrator& integrator, Scalar stepLength)
    {
        integrator.integrate(particles, stepLength);
        return stepLength;
    }
    /*void setRelativeError(Scalar relError){ relativeError = relError; }
    void setAbsoluteError(Scalar absError){ absoluteError = absError; }
private:
    Scalar absoluteError{1e-15};
    Scalar relativeError{1e-15};*/
};

#endif
