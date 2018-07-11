////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:particleSystem.h                                                                                           //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  The main class of this n-body code. This class includes                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H
#include <fstream>
#include <cstring>
#include <iomanip>

/**  @brief Base class of particle System.
 *
 *   Base particles system class. Other particle system can inherit this class. Considering the performance, we don't
 *   set virtual function.
 */
template <typename Particles, typename Interaction>
class particleSystem
{
public:
    /* Typedef */
    template<typename T, size_t S>
    using Container   = typename Particles::template Container<T, S>;
    
    using Scalar      = typename Particles::Scalar;
    
    using Vector      = typename Particles::Vector;
    
    using VectorArray = typename Particles::VectorArray;
    
    using ScalarArray = typename Particles::ScalarArray;
    
    using IntArray    = typename Particles::IntArray;
    
    using SizeArray   = typename Particles::SizeArray;
    
    using ActiveScalarArray = typename Particles::ActiveScalarArray;
    /* Typedef */
    
    constexpr static size_t arraySize{Particles::arraySize};
    
    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    inline size_t particleNumber()
    {
        return partc.particleNumber();
    }
    
    /**  @brief Physical time scalar const interface. Reference to state.time*/
    inline const Scalar& time() const { return partc.time(); }
    
    /**  @brief Position array const interface. Reference to state.pos*/
    inline const VectorArray& pos() const { return partc.pos(); }
    
    /**  @brief Velocity array const interface. Reference to state.vel*/
    inline const VectorArray& vel() const { return partc.vel(); }
    
    /**  @brief Mass array const interface. Reference to attribute.mass.*/
    inline const ScalarArray& mass() const { return partc.mass(); }
    
    /**  @brief Radius array const interface. Reference to attribute.radius.*/
    inline const ScalarArray& radius() const { return partc.radius(); }
    
    /**  @brief Particle type array const interface. Reference to attribute.type.*/
    inline const IntArray& type() const { return partc.type(); }
    
    /**  @brief Particle id array const interface. Reference to attribute.type.*/
    inline const IntArray& idn() const { return partc.idn(); }
    
    /**  @brief Position vector const interface. Reference to state.pos[i]*/
    inline const Vector& pos(size_t i) const { return partc.pos(i); }
    
    /**  @brief Velocity vecotr const interface. Reference to state.vel[i]*/
    inline const Vector& vel(size_t i) const { return partc.vel(i); }
    
    /**  @brief Mass const interface. Reference to attribute.mass[i].*/
    inline const Scalar& mass(size_t i) const { return partc.mass(i); }
    
    /**  @brief Radius const interface. Reference to attribute.radius[i].*/
    inline const Scalar& radius(size_t i) const { return partc.radius(i); }
    
    /**  @brief Particle type const interface. Reference to attribute.type[i].*/
    inline const int& type(size_t i) const { return partc.type(i); }
    
    /**  @brief Particle id const interface. Reference to attribute.type[i].*/
    inline const int& idn(size_t i) const { return partc.idn(i); }


    /** @brief Interface to rescale the time.
     *
     *  Interace used by dynamic system. Transfer integration time(For some system, integration time is different from
     *  physical time) to physical time.
     *  @return The phsyical time.
     */
    Scalar timeScale(Scalar scale)
    {
        return scale;
    }
    
    /** @brief Advance position one step with current velocity. */
    void drift(Scalar stepSize)
    {
        partc.advancePos(partc.vel(), stepSize);
        partc.advanceTime(stepSize);
    }
    
    /** @brief Advance velocity one step with current acceleration. */
    void kick(Scalar stepSize)
    {
        act.zeroTotalAcc();
        
        act.calcuVelIndepAcc(partc.mass(), partc.pos(), partc.vel());
        act.calcuExtVelIndepAcc(partc.mass(), partc.pos(), partc.vel());
        
        advanceVels<Interaction::isVelDep>(stepSize);
    }
    
    /** @brief Preprocess before iteration*/
    void preIterProcess() {}
    
    /** @brief After process after iteration*/
    void afterIterProcess() {}
    
    /** @brief Virtualize default destructor.*/
    virtual ~particleSystem() {}
    
    /** @brief Overload operator << */
    friend std::ostream& operator<<(std::ostream& os, const particleSystem& sys)
    {
        return os << sys.partc;
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& is, particleSystem& sys)
    {
        return is >> sys.partc;
    }
    
    /** @brief Input variables with plain scalar array.*/
    friend ActiveScalarArray& operator>>(ActiveScalarArray& data, particleSystem& sys)
    {
        return data >> sys.partc;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend ActiveScalarArray& operator<<(ActiveScalarArray& data, const particleSystem& sys)
    {
        return data << sys.partc;
    }
    
protected:
    /**  @brief Particle class*/
    Particles partc;
    
    /**  @brief Interaction class*/
    Interaction act;
    
private:
    template<bool isVelDep>
    inline typename std::enable_if<isVelDep==false>::type
    advanceVels(Scalar stepSize)
    {
        act.calcuTotalAcc();
        partc.advanceVel(act.totalAcc(), stepSize);
    }
    
    template<bool isVelDep>
    inline typename std::enable_if<isVelDep==true>::type
    advanceVels(Scalar stepSize)
    {
        act.calcuVelDepAcc(partc.mass(), partc.pos(), partc.vel());
        act.calcuExtVelDepAcc(partc.mass(), partc.pos(), partc.vel());
        
        act.calcuTotalAcc();
        partc.advanceAuxiVel(act.totalAcc(), stepSize*0.5);
        
        act.calcuVelDepAcc(partc.mass(), partc.pos(), partc.auxiVel());
        act.calcuExtVelDepAcc(partc.mass(), partc.pos(), partc.auxiVel());
        
        act.calcuTotalAcc();
        partc.advanceVel(act.totalAcc(), stepSize);
        
        act.calcuVelDepAcc(partc.mass(), partc.pos(), partc.vel());
        act.calcuExtVelDepAcc(partc.mass(), partc.pos(), partc.vel());
        
        act.calcuTotalAcc();
        partc.advanceAuxiVel(act.totalAcc(), stepSize*0.5);
    }
};


#endif
