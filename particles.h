
#ifndef PARTICLES_H
#define PARTICLES_H
#include "protoType.h"
#include "libs.h"
#include "devTools.h"
namespace SpaceH
{
    
/**
 *  @brief Basic velocity independent particles group.
 *  @tparam Dtype Type of scalar. e.g., float, double, kahanNumber...
 *  @tparam ArraySize The size of the arrays in whole system. SpaceH::DYNAMICAL for dynamical array.
 */
template<typename Dtype, size_t ArraySize>
class VelIndepParticles
{
public:
    /* Typedef */
    using type         = SpaceH::ProtoType<Dtype, ArraySize>;
    using Scalar       = typename type::Scalar;
    using Vector       = typename type::Vector;
    using VectorArray  = typename type::VectorArray;
    using ScalarArray  = typename type::ScalarArray;
    using IntArray     = typename type::IntArray;
    using ScalarBuffer = typename type::ScalarBuffer;
    /* Typedef */
    
    /*Template parameter check*/
    CHECK_POD(Dtype)
    /*Template parameter check*/
    
    constexpr static SpaceH::DATASTRUCT dataStruct{SpaceH::DATASTRUCT::PLAIN};
    
    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    inline size_t particleNumber() const
    {
        return pos_.size();
    }
    
    /**  @brief Physical time scalar const interface. Reference to time_*/
    inline const Scalar& time() const { return time_; }
    
    /**  @brief Position array const interface. Reference to pos_*/
    inline const VectorArray& pos() const { return pos_; }
    
    /**  @brief Velocity array const interface. Reference to vel_*/
    inline const VectorArray& vel() const { return vel_; }
    
    /**  @brief Mass array const interface. Reference to mass_*/
    inline const ScalarArray& mass() const { return mass_; }
    
    /**  @brief Radius array const interface. Reference to radius_*/
    inline const ScalarArray& radius() const { return radius_; }
    
    /**  @brief Particle type array const interface. Reference to type_.*/
    inline const IntArray& kind() const { return type_; }
    
    /**  @brief Particle id array const interface. Reference to type_.*/
    inline const IntArray& idn() const { return idn_; }
    
    /**  @brief Position vector const interface. Reference to pos_[i]*/
    inline const Vector& pos(size_t i) const { return pos_[i]; }
    
    /**  @brief Velocity vecotr const interface. Reference to vel_[i]*/
    inline const Vector& vel(size_t i) const { return vel_[i]; }
    
    /**  @brief Mass const interface. Reference to mass_[i].*/
    inline const Scalar& mass(size_t i) const { return mass_[i]; }
    
    /**  @brief Radius const interface. Reference to radius_[i].*/
    inline const Scalar& radius(size_t i) const { return radius_[i]; }
    
    /**  @brief Particle type const interface. Reference to type_[i].*/
    inline const int& kind(size_t i) const { return type_[i]; }
    
    /**  @brief Particle id const interface. Reference to type_[i].*/
    inline const int& idn(size_t i) const { return idn_[i]; }
    
    /** @brief Advance the time.
     *  @param dt Time increament.
     */
    inline void advanceTime(Scalar dt)
    {
        SpaceH::advanceScalar(time_, dt);
    }
    
    /** @brief Advance the position array with internal velocity array.
     *  @param stepSize The advance step size.
     */
    inline void advancePos(Scalar stepSize)
    {
        SpaceH::advanceVector(pos_, vel_, stepSize);
    }
    
    /** @brief Advance the  velocity array with given acceleration array.
     *  @param stepSize The advance step size.
     *  @param acc      The acceleration array.
     */
    inline void advanceVel(const VectorArray& acc, Scalar stepSize)
    {
        SpaceH::advanceVector(vel_, acc, stepSize);
    }

    /** @brief Input(Initialize) variables with istream.*/
    friend std::istream& operator>>(std::istream& is, VelIndepParticles& partc)
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
        
        Vector CMPos = SpaceH::calcuCMCoord(partc.mass_, partc.pos_, partc.totalMass_);
        Vector CMVel = SpaceH::calcuCMCoord(partc.mass_, partc.vel_, partc.totalMass_);
        
        SpaceH::moveToCMCoord(partc.pos_, CMPos);
        SpaceH::moveToCMCoord(partc.vel_, CMVel);
        
        return is;
    }
    
    /** @brief Output variables to ostream.*/
    friend std::ostream& operator<<(std::ostream& os, const VelIndepParticles& partc)
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
    friend size_t operator>>(const ScalarBuffer& data, VelIndepParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        size_t loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            partc.pos_[i].x = data[loc++];
            partc.pos_[i].y = data[loc++];
            partc.pos_[i].z = data[loc++];
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            partc.vel_[i].x = data[loc++];
            partc.vel_[i].y = data[loc++];
            partc.vel_[i].z = data[loc++];
        }
        
        partc.time_ = data[loc++];
        
        return loc;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend size_t operator<<(ScalarBuffer& data, const VelIndepParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        data.clear();
        data.reserve(particleNum*6 + 1);
        
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            data.emplace_back(partc.pos_[i].x);
            data.emplace_back(partc.pos_[i].y);
            data.emplace_back(partc.pos_[i].z);
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            data.emplace_back(partc.vel_[i].x);
            data.emplace_back(partc.vel_[i].y);
            data.emplace_back(partc.vel_[i].z);
        }
        
        data.emplace_back(partc.time_);

        return data.size();
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

/**
 *  @brief Basic velocity dependent particles group.
 *  @tparam Dtype Type of scalar. e.g., float, double, kahanNumber...
 *  @tparam ArraySize The size of the arrays in whole system. SpaceH::DYNAMICAL for dynamical array.
 */
template<typename Dtype, size_t ArraySize>
class VelDepParticles : public VelIndepParticles<Dtype, ArraySize>
{
public:
    /* Typedef */
    using Base = VelIndepParticles<Dtype, ArraySize>;
    
    using typename Base::type;
    
    using Scalar = typename type::Scalar;
    
    using Vector = typename type::Vector;
    
    using VectorArray = typename type::VectorArray;
    
    using ScalarBuffer = typename type::ScalarBuffer;
    /* Typedef */
    
    /*Template parameter check*/
    CHECK_POD(Dtype)
    /*Template parameter check*/
    
    /**  @brief Auxiliary velocity array const interface. Reference to auxi_vel_*/
    inline const VectorArray& auxiVel() const { return auxi_vel_; }
    
    /**  @brief Auxiliary velocity vecotr const interface. Reference to auxi_vel_[i] */
    inline const Vector& auxiVel(size_t i) const { return auxi_vel_[i]; }
    
    /** @brief Advance the auxiliary velocity array with given acceleration array.
     *  @param stepSize The advance step size.
     *  @param acc      The acceleration array.
     */
    inline void advanceAuxiVel(const VectorArray& acc, Scalar stepSize)
    {
        advanceVector(auxi_vel_, acc, stepSize);
    }
    
    /**  @brief synchronize auxiVel_ with vel_ */
    inline void synAuxiVelwithVel()
    {
        auxi_vel_ = this->vel_;
    }
    
    /** @brief Input(Initialize) variables with istream.*/
    friend std::istream& operator>>(std::istream& is, VelDepParticles& partc)
    {
        is >> static_cast<Base&>(partc);
        
        partc.auxi_vel_ = partc.vel_;
        return is;
    }
    
    /** @brief Input variables with plain scalar array.*/
    friend size_t operator>>(const ScalarBuffer& data, VelDepParticles& partc)
    {
        size_t loc = data >> static_cast<Base&>(partc);
        
        size_t particleNum = partc.particleNumber();
    
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            partc.auxi_vel_[i].x = data[loc++];
            partc.auxi_vel_[i].y = data[loc++];
            partc.auxi_vel_[i].z = data[loc++];
        }
        
        return loc;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend size_t operator<<(ScalarBuffer& data, const VelDepParticles& partc)
    {
        size_t loc = data << static_cast<const Base&>(partc);
        
        size_t particleNum = partc.particleNumber();
        
        data.reserve(loc + particleNum*3 );
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            data.emplace_back(partc.auxi_vel_[i].x);
            data.emplace_back(partc.auxi_vel_[i].y);
            data.emplace_back(partc.auxi_vel_[i].z);
        }
        
        return data.size();
    }
protected:
    /** @brief Auxiliary velocity array of the particles. Element is 3D vector.*/
    VectorArray auxi_vel_;
};
    
    /**
     *  @brief Basic particles group, wrapper on VelIndepParticles and VelDepParticles.
     *  @tparam Dtype Type of scalar. e.g., float, double, kahanNumber...
     *  @tparam ArraySize The size of the arrays in whole system. SpaceH::DYNAMICAL for dynamical array.
     *  @tparam IsVelDep Template parameters to determine if the particles are velocity dependent.
     */
    template<typename Dtype, size_t ArraySize, bool IsVelDep>
    struct Particles : public VelIndepParticles<Dtype, ArraySize>
    {
        constexpr static bool isVelDep{false};
    };
    
    template<typename Dtype, size_t ArraySize>
    struct Particles<Dtype, ArraySize, true> : public VelDepParticles<Dtype, ArraySize>
    {
        constexpr static bool isVelDep{true};
    };
}
#endif

