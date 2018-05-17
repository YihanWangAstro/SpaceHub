/**
 *                                          _ooOoo_
 *                                         o8888888o
 *                                         88" . "88
 *                                         (| -_- |)
 *                                          O\ = /O
 *                                      ____/`---'\____
 *                                    .   ' \\| |// `.
 *                                     / \\||| : |||// \
 *                                  / _||||| -:- |||||- \
 *                                     | | \\\ - /// | |
 *                                  | \_| ''\---/'' | |
 *                                    \ .-\__ `-` ___/-. /
 *                                ___`. .' /--.--\ `. . __
 *                              ."" '< `.___\_<|>_/___.' >'"".
 *                             | | : `- \`.;`\ _ /`;.`/ - ` : | |
 *                              \ \ `-. \_ __\ /__ _/ .-` / /
 *                      ======`-.____`-.___\_____/___.-`____.-'======
 *                                          `=---='
 *
 *                      .............................................
 *                                佛祖保佑             永无BUG
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:dynamicSystem.h                                                                                            //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DYNAMICSYSTEM_H
#define DYNAMICSYSTEM_H
#include "errhand.h"
#include <float.h>
#include <fstream>
template<typename ParticSys, typename Integrator, typename ODEiterator>
class dynamicSystem
{
    //////////////////////////////////size_terface/////////////////////////////////////
public:
    void advanceOneStep();
    void loadText(char const* initFilePath);
    void setStepLength(double);
    virtual ~dynamicSystem(){}
    ///////////////////////////////Member variables/////////////////////////////////
public:
    double      stepLength{0.0};
    ParticSys   particles;
    Integrator  integrator;
    ODEiterator iterator;
    //////////////////////////////Private Function//////////////////////////////////
private:
    void getInitStepLength();
};
    /////////////////////////////Implement Function/////////////////////////////////
template<typename ParticSys, typename Integrator, typename ODEiterator>
inline void dynamicSystem<ParticSys, Integrator, ODEiterator>::advanceOneStep()
{
    stepLength = iterator.iterate(particles, integrator, stepLength);
}

template<typename ParticSys, typename Integrator, typename ODEiterator>
void dynamicSystem<ParticSys, Integrator, ODEiterator>::getInitStepLength()
{
    if(stepLength == 0.0)
    {
        stepLength = 0.1*YEAR*particles.timeScale(1);
    }
}

template<typename ParticSys, typename Integrator, typename ODEiterator>
void dynamicSystem<ParticSys, Integrator, ODEiterator>::loadText(char const* initFilePath)
{
    std::ifstream inFile(initFilePath);
    if(inFile.is_open())
    {
        char   head;
        size_t num;
        inFile >> head >> num;
        if(head == '#' && num == particles.size())
            inFile >> particles;
        else
            throw errhand("Particle number dismatch or wrong initial file header format!", __FILE__, __LINE__);
    }
    else
        throw errhand("Fail to open the initial condition file!", __FILE__, __LINE__);
    inFile.close();
    getInitStepLength();
}

template<typename ParticSys, typename Integrator, typename ODEiterator>
void dynamicSystem<ParticSys, Integrator, ODEiterator>::setStepLength(double stepSize)
{
    stepLength = stepSize;
}

template<typename ParticSys,
         template<typename> class Integrator,
         template<typename,typename> class ODEiterator>
using spaceX = dynamicSystem<ParticSys, Integrator<ParticSys>, ODEiterator<ParticSys, Integrator<ParticSys> > >;
#endif

