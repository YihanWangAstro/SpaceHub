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
 */
#ifndef DYNAMICSYSTEM_H
#define DYNAMICSYSTEM_H

#include "devTools.h"
#include <fstream>
#include <vector>
#include "macros.h"

namespace SpaceH {
/**
 *
 * @tparam ParticleSys
 */
    template<typename ParticleSys>
    class RunArgs {
    public:
        using Scalar = typename ParticleSys::Scalar;
        using Callback = std::function<void (ParticleSys &)>;

        Scalar endTime{0};

        void preOption(ParticleSys &partc_sys) const {
            for (const auto &opt : preOpts) {
                opt(partc_sys);
            }
        }

        void postOption(ParticleSys &partc_sys) const {
            for (const auto &opt : postOpts) {
                opt(partc_sys);
            }
        }

        void registerPreOption(Callback&& fun) {
            preOpts.emplace_back(fun);
        }

        void registerPostOption(Callback&& fun) {
            postOpts.emplace_back(fun);
        }

    private:
        std::vector<Callback> preOpts;
        std::vector<Callback> postOpts;
    };

/**
 *  @brief A wrapper to make particle system, integrator and ODE iterator work together.
*/
    template<typename ParticSys, typename ODEiterator>
    class Solver {
    public:
        /* Typedef */
        using type = typename ParticSys::type;
        using Scalar = typename type::Scalar;
        using ScalarBuffer = typename type::ScalarBuffer;
        using RunArgs = RunArgs<ParticSys>;
        /* Typedef */

        void advanceOneStep();

        void loadText(char const *initFilePath);

        /**  @brief Set the step length*/
        inline void setStepLength(Scalar step_size) {
            stepLength = step_size;
        }

        inline Scalar getCurrentStepLength() {
            return stepLength;
        }

        void run(const RunArgs &arg) {
            Scalar end_time = arg.endTime;
            DEBUG_MSG(true, end_time);
            for (; particles.time() < end_time;) {
                arg.preOption(particles);
                advanceOneStep();
                arg.postOption(particles);
            }
        }

        virtual ~Solver() {} /**< @brief Default destructor, virtualize for inherent class*/
    public:
        /** @brief Macro step size for ODE iterator*/
        Scalar stepLength{0.0};

        /** @brief Steps*/
        int steps{0};

        /** @brief Particle system*/
        ParticSys particles;

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
    template<typename ParticSys, typename ODEiterator>
    inline void Solver<ParticSys, ODEiterator>::advanceOneStep() {
        particles.preIterProcess();
        stepLength = iterator.iterate(particles, stepLength);
        particles.afterIterProcess();
        steps++;
    }

/**  @brief Calculate the initial step length of the particle system
 *
 *   If the user didn't set the step length with setStepLength(), calculate the proper initial
 *   step length automatically.
 *
 */
    template<typename ParticSys, typename ODEiterator>
    void Solver<ParticSys, ODEiterator>::getInitStepLength() {
        if (stepLength == 0.0) {
            stepLength = 0.1 * particles.timeScale();
            //std::cout << "dyn T=" << stepLength/YEAR << '\n';
        }
        steps = 0;
        particles.evaluateAcc();
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
    template<typename ParticSys, typename ODEiterator>
    void Solver<ParticSys, ODEiterator>::loadText(char const *initFilePath) {
        std::ifstream inFile(initFilePath);

        if (inFile.is_open()) {
            char header;
            inFile >> header;
            if (header == '#') {
                inFile >> particles;
                inFile.close();
            } else {
                inFile.close();
                SpaceH::errMsg("Input file header should begin with '#'.", __FILE__, __LINE__);
            }
        } else {
            inFile.close();
            SpaceH::errMsg("Fail to open the initial condition file!", __FILE__, __LINE__);
        }

        getInitStepLength();
    }

}
#endif

