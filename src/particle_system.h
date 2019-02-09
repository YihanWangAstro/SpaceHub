
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "core_computation.h"
#include "dev_tools.h"
#include "macros.h"
#include "type_class.h"

namespace SpaceH {


    template<typename Particles, typename Acc, typename EoM>
    class ParticleSystem {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        inline size_t number() {
            return partc_.number();
        }

        void advanceTime(Scalar dt) {
            partc_.time_ += dt;
        }

        void advancePos(Scalar stepSize) {
            advancePos(partc_.vx, partc_.vy, partc_.vz, stepSize);
        }

        void advancePos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            comp::advanceVector(partc_.px, partc_.vx, stepSize);
            comp::advanceVector(partc_.py, partc_.vy, stepSize);
            comp::advanceVector(partc_.pz, partc_.vz, stepSize);
        }

        void advanceVel(Scalar stepSize) {
            advanceVel(acc_.ax, acc_.ay, acc_.az, stepSize);
        }

        void advanceVel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            comp::advanceVector(partc_.vx, ax, stepSize);
            comp::advanceVector(partc_.vy, ay, stepSize);
            comp::advanceVector(partc_.vz, az, stepSize);
        }

        void evaluateAcc() {
            eom_.eval_acc(partc_, acc_);
        }

        void drift(Scalar stepSize) {
            advanceTime(stepSize);
            advancePos(stepSize);
        }

        void kick(Scalar stepSize) {

        }
    private:
        Particles partc_;
        Acc acc_;
        EoM eom_;
    };

/**  @brief Base class of particle System.
 *
 *   Base particles system class. Other particle system can inherit this class. Considering the performance, we don't
 *   set virtual function, derived class may hide the base class function name. Therefore, be careful about downcast.
 *   Further test will be performed to virtualize non-hot functions.
 */
    template<typename Particles, typename Interaction>
    class ParticleSystem {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);
        using State        = typename Particles::State;
        using ParticleType = Particles;
        /* Typedef */

        /** @brief interface adapter to inherit the interface of the data member
         *  The macros take four args (TYPE, MEMBER, NAME, NEWNAME). Each macros create two interfaces, they are:
         *
         *  1. const TYPE &NEWNAME () const { return MEMBER.NAME();};
         *  2. const typename TYPE::value_type & NEWNAME (size_t i) const { return MEMBER.NAME(i);};
         *  See macros definition in 'devTools.h'
         */
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(pos,    VectorArray, partc, pos   );
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(vel,    VectorArray, partc, vel   );
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(mass,   ScalarArray, partc, mass  );
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(radius, ScalarArray, partc, radius);
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(idn,    IndexArray,  partc, idn   );
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_SCALAR(time,  Scalar,      partc, time  );
        SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(acc,    VectorArray, act,   acc   );

        /*Template parameter check*/
        CHECK_TYPE(Particles, Interaction);
        /*Template parameter check*/

        constexpr static bool isVelDep{Interaction::isVelDep};

        /** @brief Get the number of the particles.
         *  @return The particle number.
         */
        inline size_t particleNumber() const {
            return partc.particleNumber();
        }


        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void drift(Scalar stepSize) {
            partc.advancePos(stepSize);
            partc.advanceTime(stepSize);
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            act.zeroTotalAcc();
            act.calcuVelIndepAcc(partc);

            /*after long time struggling, decide to use 'constexpr if' in c++17 to improve readability. The price is low
            version compiler will treat this as normal 'if' such that unrelavent brach will be checked at runtime.
            Although, branch prediction somewhat alleviate this overhead, I would suspect the instructions of
            unrelevant branch that generated at compile time may affect the efficency of CPU pipline.*/
            if constexpr (!Interaction::isVelDep){
                act.sumTotalAcc();
                partc.advenceVel(act.acc(), stepSize);
            }else {
                State v0 = partc.vel_state();
                act.calcuVelDepAcc(partc);
                act.sumTotalAcc();
                partc.advenceVel(act.acc(), 0.5*stepSize);
                iterateVeltoConvergent(v0, 0.5*stepSize);
                partc.advenceVel(act.acc(), 0.5*stepSize);
            }
        }

        inline void advanceTime(Scalar dt) {
            partc.advenceTime(dt);
        }

        /** @brief Advance the position array with internal velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(Scalar stepSize) {
            partc.advancePos(stepSize);
        }

        /** @brief Advance the position array with given velocity array.
         *  @param vel The given velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(const VectorArray &vel, Scalar stepSize) {
            partc.advancePos(vel,stepSize);
        }

        /** @brief Advance the  velocity array with given acceleration array.
         *  @param stepSize The advance step size.
         *  @param acc      The acceleration array.
         */
        inline void advanceVel(const VectorArray &acc, Scalar stepSize) {
            partc.advanceVel(acc,stepSize);
        }

        /** @brief Evaluate acceleration with current velocity and position without any advance
         *  @return The acceleration after evaluation.
         *  @note Used in non-symplectic method for function evaluation
         */
        const VectorArray& evaluateAcc() {
            act.zeroTotalAcc();
            act.calcuVelIndepAcc(partc);
            act.calcuVelDepAcc(partc);
            act.sumTotalAcc();

            return act.acc();
        }

        /** @brief Preprocess before iteration*/
        void preIterProcess() {}

        /** @brief After process after iteration*/
        void afterIterProcess() {
            partc.update();
        }

        /** @brief Virtualize default destructor.*/
        virtual ~ParticleSystem() {}

        /** @brief Overload operator << */
        friend std::ostream &operator<<(std::ostream &os, const ParticleSystem &sys) {
            sys.writeHeader(os);
            sys.partc.write(os, SpaceH::Unit::STD_UNIT);
            return os;
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &is, ParticleSystem &sys) {
            size_t num = sys.readHeader(is);

            if constexpr (ParticleSystem::Types::array_size == SpaceH::DYNAMICAL) {
                sys.resize(num);
            }
            if (num == sys.particleNumber()) {
                sys.partc.read(is, SpaceH::Unit::STD_UNIT);
                sys.partc.moveToCoM();
            } else {
                SPACEHUB_ABORT("You are using fixed particle number system, the particle number in initial file is not consistent with the system you are using!");
            }
            return is;
        }

        /** @brief Input variables with plain scalar array.*/
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {
            return partc.read(data, flag);
        }

        /** @brief Output variables to plain scalar array.*/
        size_t write(ScalarBuffer &data, const IO_flag flag = IO_flag::STD) const {
            return partc.write(data, flag);
        }

        Vector posCoM() {
            return SpaceH::calcuCMCoord(partc.mass(), partc.pos());
        }
        Vector velCoM() {
            return SpaceH::calcuCMCoord(partc.mass(), partc.vel());
        }
    protected:
        /**  @brief Particle class*/
        Particles partc;

        /**  @brief Interaction class*/
        Interaction act;

    private:
        void writeHeader(std::ostream &os) const {
            os << '#' << partc.particleNumber();
        }

        size_t readHeader(std::istream &is) {
            char tag;
            is >> tag;

            if (tag == '#') {
                size_t particleNum;
                is  >> particleNum;
                return particleNum;
            } else {
                SPACEHUB_ABORT("Input file header should begin with '#'.");
            }
        }

        bool isConvergent(const VectorArray & v1, const VectorArray & v2) {
            size_t size = particleNumber();
            Scalar max_dv = 0;
            for(size_t i = 0 ; i < size; ++i) {
                Scalar dv = SpaceH::max(max_dv, ((v1[i] - v2[i])/v2[i]).abs().max_component());
            }
            return max_dv < SpaceH::epsilon<Scalar>::value;
        }

    protected:
        void iterateVeltoConvergent(const State &v0, Scalar stepSize) {
            for(size_t  i = 0 ; i < 10; ++i) {
                State v_new = partc.vel_state();
                act.calcuVelDepAcc(partc);
                act.sumTotalAcc();
                partc.set_vel_state(v0);
                partc.advenceVel(act.acc(), stepSize);
                if(isConvergent(v_new.raw(), partc.vel_state().raw()))
                    return;
            }
        }
    };
}

#endif
