/*
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
 */
#ifndef DYNAMICSYSTEM_H
#define DYNAMICSYSTEM_H
#include "errhand.h"
#include <float.h>
#include <fstream>

/**
 *  @brief A wrapper to make particle system, integrator and ODE iterator work together.
*/
template<typename ParticSys, typename Integrator, typename ODEiterator>
class dynamicSystem
{
public:
    void advanceOneStep();
    void loadText(char const* initFilePath);
    void setStepLength(double);
    virtual ~dynamicSystem() {} /**< @brief Default destructor, virtualize for inherent class*/
public:
    /** @brief Macro step size for ODE iterator*/
    double stepLength{0.0};
    
    /** @brief Particle system*/
    ParticSys particles;
    
    /** @brief Integrator*/
    Integrator integrator;
    
    /** @brief ODE Iterator*/
    ODEiterator iterator;
    
private:
    void getInitStepLength();
};

/**  @brief Advance the particle system for one step.
 *
 *   Advance the particle system with current steplength stepLength. The ODE iterator
 *   iterate the integrator to convergence by its own implement. The step length will also
 *   be updated by its own implement.
 */
template<typename ParticSys, typename Integrator, typename ODEiterator>
inline void dynamicSystem<ParticSys, Integrator, ODEiterator>::advanceOneStep()
{
    stepLength = iterator.iterate(particles, integrator, stepLength);
}

/**  @brief Calculate the initial step length of the particle system
 *
 *   If the user didn't set the step length with setStepLength(), calculate the proper initial
 *   step length automatically.
 *
 */
template<typename ParticSys, typename Integrator, typename ODEiterator>
void dynamicSystem<ParticSys, Integrator, ODEiterator>::getInitStepLength()
{
    if(stepLength == 0.0)
    {
        stepLength = 0.1 * YEAR * particles.timeScale(1);
    }
}


/**  @brief Load particle system initial condition from file
 *
 *   This function will read and check the initial file header (begin with '#') and the
 *   particle number after the '#'. Pass the rest information to particles by
 *   operator '>>'. The way to load the initial condition depend on the implemet of
 *   the particles. If the initial condition read successfully. This function will
 *   call getInitStepLength() to set the initial step length.
 *
 *   @param initFilePath The relative path of initial conditions file
 *   @exception If the partcile number in the header is inconsisitent with the size of
 *              particles, this function will throw an exception.
 */
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

/**  @brief Set the step length*/
template<typename ParticSys, typename Integrator, typename ODEiterator>
void dynamicSystem<ParticSys, Integrator, ODEiterator>::setStepLength(double stepSize)
{
    stepLength = stepSize;
}

/**  @brief Alias of template name, linking the particle system, integrator and ODE iterator*/
template<typename ParticSys,
         template<typename> class Integrator,
         template<typename, typename> class ODEiterator>
using spaceX = dynamicSystem<ParticSys, Integrator<ParticSys>, ODEiterator<ParticSys, Integrator<ParticSys>>>;
#endif

