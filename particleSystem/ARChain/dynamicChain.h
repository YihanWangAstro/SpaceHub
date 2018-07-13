
#ifndef DYNAMICCHAIN_H
#define DYNAMICCHAIN_H
#include "../../particles.h"
#include "chain.h"
/**
 *  @brief Class of dynamical variable.
 *
 */
template<typename Dtype, size_t ArraySize>
class ChainParticles : public ReguParticles<Dtype, ArraySize>
{
public:
    /* Typedef */
    using Base = ReguParticles<Dtype, ArraySize>;
    
    using typename Base::type;
    
    using Scalar = typename type::Scalar;
    
    using Vector = typename type::Vector;
    
    using VectorArray = typename type::VectorArray;
    
    using IndexArray = typename type::IndexArray;
    
    using ActiveScalarArray = typename Base::ActiveScalarArray;
    /* Typedef */
    
    constexpr static SpaceH::DATASTRUCT dataStruct{SpaceH::DATASTRUCT::CHAIN};
    
    using Base::mass_;
    using Base::pos_;
    using Base::vel_;
    
    /**  @brief Position array const interface. Reference to pos_*/
    inline const VectorArray& chainPos() const { return chain_pos_; }
    
    /**  @brief Velocity array const interface. Reference to vel_*/
    inline const VectorArray& chainVel() const { return chain_vel_; }
    
    /**  @brief Index array const interface. Reference to ch_index_*/
    inline const IndexArray& chainIndex() const { return ch_index_; }
    
    /**  @brief Position vector const interface. Reference to pos_[i]*/
    inline const Vector& chainPos(size_t i) const { return chain_pos_[i]; }
    
    /**  @brief Velocity vecotr const interface. Reference to vel_[i]*/
    inline const Vector& chainVel(size_t i) const { return chain_vel_[i]; }
    
    /**  @brief Index const interface. Reference to ch_index_[i]*/
    inline const size_t chainIndex(size_t i) const { return ch_index_[i]; }
    
    
    /** @brief Advance the position array with internal velocity array.
     *  @param stepSize The advance step size.
     */
    inline void advancePos(Scalar stepSize)
    {
        SpaceH::advanceVector(chain_pos_, chain_vel_, stepSize);
        SpaceH::chain::synCartesian(chain_pos_, pos_, ch_index_);
        SpaceH::moveToCMCoord(mass_, pos_);
    }
    
    /** @brief Advance the  velocity array with given acceleration array.
     *  @param stepSize The advance step size.
     *  @param acc      The acceleration array.
     */
    inline void advanceVel(const VectorArray& acc, Scalar stepSize)
    {
        const size_t chain_num = this->particleNumber() - 1;
        
        VectorArray chain_acc = acc;
        
        for(size_t i = 0 ; i < chain_num; ++i)
            chain_acc[i] = acc[ch_index_[i + 1]] - acc[ch_index_[i]];
        
        SpaceH::advanceVector(chain_vel_, chain_acc, stepSize);
        SpaceH::chain::synCartesian(chain_vel_, vel_, ch_index_);
        SpaceH::moveToCMCoord(mass_, vel_);
    }
    
    /** @brief Input(Initialize) variables with istream.*/
    friend std::istream& operator>>(std::istream& is, ChainParticles& partc)
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
        
        SpaceH::moveToCMCoord(partc.mass_,partc.pos_);
        SpaceH::moveToCMCoord(partc.mass_,partc.vel_);
        
        partc.omega_ = partc.getCapitalOmega();
        partc.bindE_ = -getTotalEnergy(partc.mass(), partc.pos(), partc.vel());
        
        SpaceH::chain::getChainIndex(partc.pos_,  partc.ch_index_);
        SpaceH::chain::synChain(partc.pos_, partc.chain_pos_, partc.ch_index_);
        SpaceH::chain::synChain(partc.vel_, partc.chain_vel_, partc.ch_index_);
        
        return is;
    }
    
    /** @brief Input variables with plain scalar array.*/
    friend ActiveScalarArray& operator>>(ActiveScalarArray& data, ChainParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        size_t d_loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            partc.chain_pos_[i].x = data[d_loc++];
            partc.chain_pos_[i].y = data[d_loc++];
            partc.chain_pos_[i].z = data[d_loc++];
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            partc.chain_vel_[i].x = data[d_loc++];
            partc.chain_vel_[i].y = data[d_loc++];
            partc.chain_vel_[i].z = data[d_loc++];
        }
        
        partc.time_  = data[d_loc++];
        partc.omega_ = data[d_loc++];
        partc.bindE_ = data[d_loc];
        
        SpaceH::chain::synCartesian(partc.chain_vel_, partc.vel_, partc.ch_index_);
        SpaceH::chain::synCartesian(partc.chain_pos_, partc.pos_, partc.ch_index_);
        SpaceH::moveToCMCoord(partc.mass_,partc.pos_);
        SpaceH::moveToCMCoord(partc.mass_,partc.vel_);
        
        return data;
    }
    
    /** @brief Output variables to plain scalar array.*/
    friend ActiveScalarArray& operator<<(ActiveScalarArray& data, const ChainParticles& partc)
    {
        size_t particleNum = partc.particleNumber();
        
        size_t d_loc = 0;
        //for locality, split into two loops
        for(size_t i = 0; i < particleNum; ++i)
        {
            data[d_loc++] = partc.chain_pos_[i].x;
            data[d_loc++] = partc.chain_pos_[i].y;
            data[d_loc++] = partc.chain_pos_[i].z;
        }
        
        for(size_t i = 0 ; i < particleNum; ++i)
        {
            data[d_loc++] = partc.chain_vel_[i].x;
            data[d_loc++] = partc.chain_vel_[i].y;
            data[d_loc++] = partc.chain_vel_[i].z;
        }
        
        data[d_loc++] = partc.time_;
        data[d_loc++] = partc.omega_;
        data[d_loc]   = partc.bindE_;
        
        return data;
    }

private:
    /** @brief Chained position array of the particles. Element is 3D vector.*/
    VectorArray chain_pos_;
    
    /** @brief Chained velocity array of the particles. Element is 3D vector.*/
    VectorArray chain_vel_;
    
    /** @brief Index array from Cartesian coordinates to chain coordinates Element is 3D vector.*/
    IndexArray ch_index_;
};
#endif

