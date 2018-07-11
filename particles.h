
#ifndef PARTICLES_H
#define PARTICLES_H
#include "protoType.h"
#include "libs.h"
/**
 *  @brief Class of dynamical variable.
 *
 */
template<typename Dtype, size_t ArraySize>
class particles : public type::ProtoType<Dtype, ArraySize>
{
public:
    /* Typedef */
    using Base = type::ProtoType<Dtype, ArraySize>;
    
    template<typename T, size_t S>
    using Container   = typename Base::template Container<T, S>;
    
    using Scalar      = typename Base::Scalar;
    
    using Vector      = typename Base::Vector;
    
    using VectorArray = typename Base::VectorArray;
    
    using ScalarArray = typename Base::ScalarArray;
    
    using IntArray    = typename Base::IntArray;
    
    using SizeArray   = typename Base::SizeArray;
    /* Typedef */
    
    constexpr static size_t activeScalar{6*Base::arraySize + 1};
    
    using ActiveScalarArray = Container<Scalar, activeScalar>;
    
    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    inline size_t particleNumber() const
    {
        return pos_.size();
    }
    
    /**  @brief Physical time scalar const interface. Reference to state.time*/
    inline const Scalar& time() const { return time_; }
    
    /**  @brief Position array const interface. Reference to state.pos*/
    inline const VectorArray& pos() const { return pos_; }
    
    /**  @brief Velocity array const interface. Reference to state.vel*/
    inline const VectorArray& vel() const { return vel_; }
    
    /**  @brief Mass array const interface. Reference to attribute.mass.*/
    inline const ScalarArray& mass() const { return mass_; }
    
    /**  @brief Radius array const interface. Reference to attribute.radius.*/
    inline const ScalarArray& radius() const { return radius_; }
    
    /**  @brief Particle type array const interface. Reference to attribute.type.*/
    inline const IntArray& type() const { return type_; }
    
    /**  @brief Particle id array const interface. Reference to attribute.type.*/
    inline const IntArray& idn() const { return idn_; }
    
    /**  @brief Position vector const interface. Reference to state.pos[i]*/
    inline const Vector& pos(size_t i) const { return pos_[i]; }
    
    /**  @brief Velocity vecotr const interface. Reference to state.vel[i]*/
    inline const Vector& vel(size_t i) const { return vel_[i]; }
    
    /**  @brief Mass const interface. Reference to attribute.mass[i].*/
    inline const Scalar& mass(size_t i) const { return mass_[i]; }
    
    /**  @brief Radius const interface. Reference to attribute.radius[i].*/
    inline const Scalar& radius(size_t i) const { return radius_[i]; }
    
    /**  @brief Particle type const interface. Reference to attribute.type[i].*/
    inline const int& type(size_t i) const { return type_[i]; }
    
    /**  @brief Particle id const interface. Reference to attribute.type[i].*/
    inline const int& idn(size_t i) const { return idn_[i]; }
    
    /** @brief Advance the time.
     *  @param dt Time increament.
     */
    inline void advanceTime(Scalar dt)
    {
        advanceScalar(time_, dt);
    }
    
    /** @brief Advance the position array with internal velocity array.
     *  @param stepSize The advance step size.
     *  @param vel      The velocity array to advance position array.
     */
    inline void advancePos(const VectorArray& vel, Scalar stepSize)
    {
        advanceVector(pos_, vel, stepSize);
    }
    
    /** @brief Advance the  velocity array with given acceleration array.
     *  @param stepSize The advance step size.
     *  @param acc      The acceleration array.
     */
    inline void advanceVel(const VectorArray& acc, Scalar stepSize)
    {
        advanceVector(vel_, acc, stepSize);
    }

    /** @brief Input(Initialize) variables with istream.*/
    friend std::istream& operator>>(std::istream& is, particles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        is >> partc.time_;
        
        partc.totalMass_ = 0;
        
        for(size_t i = 0 ; i < particleNum ; ++i)
        {
            is >> partc.idn_[i]
               >> partc.type_[i]
               >> partc.mass_[i]
               >> partc.radius_[i]
               >> partc.pos_[i]
               >> partc.vel_[i];
            
            partc.totalMass_ += partc.mass_[i];
        }
        return is;
    }
    
    /** @brief Output variables to ostream.*/
    friend std::ostream& operator<<(std::ostream& os, const particles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        os << "#" << particleNum << " " << partc.time_ << "\r\n";
        
        for(size_t i = 0 ; i < particleNum ; ++i)
        {
            os << partc.idn_[i]    << " "
               << partc.type_[i]   << " "
               << partc.mass_[i]   << " "
               << partc.radius_[i] << " "
               << partc.pos_[i]    << " "
               << partc.vel_[i]    << "\r\n";
        }
        
        return os;
    }
    
    /** @brief Input variables with plain scalar array.*/
    friend ActiveScalarArray& operator>>(ActiveScalarArray& data, particles& partc)
    {
        size_t particleNum = partc.particleNumber();
        size_t d_loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            partc.pos_[i].x = data[d_loc++];
            partc.pos_[i].y = data[d_loc++];
            partc.pos_[i].z = data[d_loc++];
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            partc.vel_[i].x = data[d_loc++];
            partc.vel_[i].y = data[d_loc++];
            partc.vel_[i].z = data[d_loc++];
        }
        
        partc.time_ = data[d_loc];
        
        return data;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend ActiveScalarArray& operator<<(ActiveScalarArray& data, const particles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        size_t d_loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            data[d_loc++] = partc.pos_[i].x;
            data[d_loc++] = partc.pos_[i].y;
            data[d_loc++] = partc.pos_[i].z;
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            data[d_loc++] = partc.vel_[i].x;
            data[d_loc++] = partc.vel_[i].y;
            data[d_loc++] = partc.vel_[i].z;
        }
        
        data[d_loc] = partc.time_;

        return data;
    }
    
protected:
    /** @brief Position array of the particles. Element is 3D vector.*/
    VectorArray pos_;
    
    /** @brief Velocity array of the particles. Element is 3D vector.*/
    VectorArray vel_;
    
    /** @brief Mass array of the particles. Element is Scalar.*/
    ScalarArray mass_;
    
    /** @brief Radius array of the particles. Element is Scalar.*/
    ScalarArray radius_;
    
    /** @brief Type Array of the particles. Element is int.*/
    IntArray type_;
    
    /** @brief Id Array of the particles. Element is int.*/
    IntArray idn_;
    
    /** @brief The physical time of the dynamic system*/
    Scalar time_;
    
    /** @brief The total mass of the system*/
    Scalar totalMass_;
};
#endif

