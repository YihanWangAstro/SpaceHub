
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "coreComputation.h"
#include "devTools.h"

namespace SpaceH {
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
        using type         = typename Particles::type;
        using Scalar       = typename type::Scalar;
        using Vector       = typename type::Vector;
        using VectorArray  = typename type::VectorArray;
        using ScalarArray  = typename type::ScalarArray;
        using IntArray     = typename type::IntArray;
        using ScalarBuffer = typename type::ScalarBuffer;
        using ParticleType = Particles;
        /* Typedef */

        /*Template parameter check*/
        CHECK_TYPE(Particles, Interaction);
        /*Template parameter check*/

        constexpr static bool isVelDep{Interaction::isVelDep};
        constexpr static size_t arraySize{type::arraySize};

        /** @brief Get the number of the particles.
         *  @return The particle number.
         */
        inline size_t particleNumber() const {
            return partc.particleNumber();
        }

        /** @brief Resize all containers if they are dynamical
         *  @param New size of container.
         */
        void resize(size_t new_siz) {
            partc.resize(new_siz);
            act.resize(new_siz);
        }

        /** @brief Reserve space for all containers if they are dynamical
         *  @param New capacity of container.
         */
        void reserve(size_t new_cap) {
            partc.reserve(new_cap);
            act.reserve(new_cap);
        }

        /**  @brief Physical time scalar const interface. Reference to partc.time*/
        inline const Scalar &time() const {
            return partc.time;
        }

        /**  @brief Position array const interface. Reference to partc.pos*/
        inline const VectorArray &pos() const {
            return partc.pos;
        }

        /**  @brief Velocity array const interface. Reference to partc.vel*/
        inline const VectorArray &vel() const {
            return partc.vel;
        }

        /**  @brief Acceleration array const interface. Reference to partc.pos*/
        inline const VectorArray &acc() const {
            return act.totalAcc();
        }

        /**  @brief Mass array const interface. Reference to partc.mass.*/
        inline const ScalarArray &mass() const {
            return partc.mass;
        }

        /**  @brief Radius array const interface. Reference to partc.radius.*/
        inline const ScalarArray &radius() const {
            return partc.radius;
        }

        /**  @brief Particle id array const interface. Reference to partc.type.*/
        inline const IntArray &idn() const {
            return partc.idn;
        }

        /**  @brief Position vector const interface. Reference to partc.pos[i]*/
        inline const Vector &pos(size_t i) const {
            return partc.pos[i];
        }

        /**  @brief Velocity vecotr const interface. Reference to partc.vel[i]*/
        inline const Vector &vel(size_t i) const {
            return partc.vel[i];
        }

        /**  @brief Acceleration vecotr const interface. Reference to partc.vel[i]*/
        inline const Vector &acc(size_t i) const {
            return act.totalAcc(i);
        }

        /**  @brief Mass const interface. Reference to partc.mass[i].*/
        inline const Scalar &mass(size_t i) const {
            return partc.mass[i];
        }

        /**  @brief Radius const interface. Reference to partc.radius[i].*/
        inline const Scalar &radius(size_t i) const {
            return partc.radius[i];
        }

        /**  @brief Particle id const interface. Reference to partc.type[i].*/
        inline const int &idn(size_t i) const {
            return partc.idn[i];
        }

        /** @brief Interface to rescale the time.
         *
         *  Interace used by dynamic system. Transfer integration time(For some system, integration time is different from
         *  physical time) to physical time.
         *  @return The phsyical time.
         */
        Scalar timeScale() {
            if (isAllZero(partc.vel))
                return SpaceH::minfallFreeTime(partc.mass, partc.pos);
            else
                return SpaceH::minAccdot(partc.mass, partc.pos, partc.vel);
        }

        /** @brief Advance position one step with current velocity. Used for symplectic integrator.*/
        void drift(Scalar stepSize) {
            SpaceH::advanceVector(partc.pos, partc.vel, stepSize);
            partc.Time += stepSize;
        }

        /** @brief Advance velocity one step with current acceleration. Used for symplectic integrator.*/
        void kick(Scalar stepSize) {
            act.zeroTotalAcc();

            act.calcuVelIndepAcc(partc);
            act.calcuExtVelIndepAcc(partc);

            /*after long time struggling, decide to use 'constexpr if' in c++17 to improve readability. The price is low
            version compiler will treat this as normal 'if' such that unrelavent brach will be checked at runtime.
            Although, branch prediction somewhat alleviate this overhead, I would suspect the instructions of
            unrelevant branch that generated at compile time may affect the efficency of CPU pipline.*/
            if constexpr (!Interaction::isVelDep){
                act.calcuTotalAcc();
                SpaceH::advanceVector(partc.vel, act.totalAcc(), stepSize);
            }else {
                Vector v0 = partc.vel;
                act.calcuVelDepAcc(partc);
                act.calcuExtVelDepAcc(partc);
                act.sumTotalAcc();
                SpaceH::advanceVector(partc.vel, act.totalAcc(), 0.5*stepSize);
                for(size_t  i = 0 ; i < Max_Iter; ++i) {
                    Vector v_half = partc.vel;
                    act.calcuVelDepAcc(partc);
                    act.calcuExtVelDepAcc(partc);
                    act.sumTotalAcc();
                    partc.vel = v0;
                    SpaceH::advanceVector(partc.vel, act.totalAcc(), 0.5*stepSize);
                    if(isConvergent(v_half, partc.vel))
                        break;
                }
                SpaceH::advanceVector(partc.vel, act.totalAcc(), 0.5*stepSize);
            }
        }

        inline void advanceTime(Scalar dt) {
            partc.time += dt;
        }

        /** @brief Advance the position array with internal velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(Scalar stepSize) {
            SpaceH::advanceVector(partc.pos, partc.vel, stepSize);
        }

        /** @brief Advance the position array with given velocity array.
         *  @param vel The given velocity array.
         *  @param stepSize The advance step size.
         */
        inline void advancePos(const VectorArray &vel, Scalar stepSize) {
            SpaceH::advanceVector(partc.pos, vel, stepSize);
        }

        /** @brief Advance the  velocity array with given acceleration array.
         *  @param stepSize The advance step size.
         *  @param acc      The acceleration array.
         */
        inline void advanceVel(const VectorArray &acc, Scalar stepSize) {
            SpaceH::advanceVector(partc.vel, acc, stepSize);
        }

        /** @brief Evaluate acceleration with current velocity and position without any advance
         *  @return The acceleration after evaluation.
         *  @note Used in non-symplectic method for function evaluation
         */
        const VectorArray& evaluateAcc() {
            act.zeroTotalAcc();
            act.calcuVelIndepAcc(partc);
            act.calcuExtVelIndepAcc(partc);
            act.calcuVelDepAcc(partc);
            act.calcuExtVelDepAcc(partc);
            act.sumTotalAcc();

            return act.totalAcc();
        }

        /** @brief Calculate the potential energy of the system*/
        inline Scalar potentialEnergy() const {
            return SpaceH::getPotentialEnergy(partc.mass, partc.pos);
        }

        /** @brief Calculate the kinetic energy of the system*/
        inline Scalar kineticEnergy() const {
            return SpaceH::getKineticEnergy(partc.mass, partc.vel);
        }

        /** @brief Calculate the total energy of the system*/
        inline Scalar totalEnergy() const {
            return potentialEnergy() + kineticEnergy();
        }

        /** @brief Preprocess before iteration*/
        void preIterProcess() {}

        /** @brief After process after iteration*/
        void afterIterProcess() {}

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

            if constexpr (sys::arraySize == SpaceH::DYNAMICAL) {
                sys.resize(num);
            }
            if (num == sys::particleNumber()) {
                sys.partc.read(is, SpaceH::Unit::STD_UNIT);
                total_mass = SpaceH::sumArray(partc.mass);
                SpaceH::moveToCoM(partc.mass, partc.pos, total_mass);
                SpaceH::moveToCoM(partc.mass, partc.vel, total_mass);
            } else {
                SpaceH::errMsg(
                        "You are using fixed particle number system, the particle number in initial file is not consistent with the system you are using!",
                        __FILE__, __LINE__);
            }
            return is;
        }

        /** @brief Input variables with plain scalar array.*/
        size_t read(const ScalarBuffer &data, const IO_flag flag = IO_flag::STD) {
            return partc.read(data, flag);
        }

        /** @brief Output variables to plain scalar array.*/
        size_t write(ScalarBuffer &data, const IO_flag IO_flag = IO_flag::STD) const {
            return partc.write(data, flag);
        }

    protected:
        /**  @brief Particle class*/
        Particles partc;

        /**  @brief Interaction class*/
        Interaction act;

        /**  @brief total mass of particles*/
        Scalar total_mass;

    private:
        void writeHeader(std::ostream &os) {
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
                SpaceH::errMsg("Input file header should begin with '#'.", __FILE__, __LINE__);
            }
        }

        bool isConvergent(Vector& v1, Vector& v2) {
            size_t size = particleNumber();
            Scalar max_dv = 0;
            for(size_t i = 0 ; i < size; ++i) {
                Scalar dv = SpaceH::max(max_value, ((v1[i] - v2[i])/v2[i]).abs().max_component());
            }
            return max_dv < SpaceH::epsilon<Scalar>::value;
        }
    };
}

#endif
