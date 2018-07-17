////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:size_tercation.h                                                                                              //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef FORCES_H
#define FORCES_H
#include "../protoType.h"
namespace SpaceH
{
    
template<typename Dtype, size_t ArraySize>
struct EmptyForce
{
    /* Typedef */
    using type = SpaceH::ProtoType<Dtype, ArraySize>;
    
    using Scalar = typename type::Scalar;
    
    using Vector = typename type::Vector;
    
    using VectorArray = typename type::VectorArray;
    
    using ScalarArray = typename type::ScalarArray;
    
    using IndexArray = typename type::IndexArray;
    /* Typedef */
    
public:
    void calcuAcc(const ScalarArray& mass, const VectorArray& pos) {}
    void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel) {}
    void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& chainPos, const IndexArray& chainIndex) {}
    void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel,
                    const VectorArray& chainPos, const VectorArray& chainVel, const IndexArray& chainIndex) {}
    
    void addTotal(VectorArray& acc) {}
    
    const VectorArray& acc() {std::cout << "empty force unaccessible data!" << "\r\n"; }
    
    const Vector& acc(size_t i) {std::cout << "empty force unaccessible data!" << "\r\n"; }
    
    constexpr static bool isVelDep{false};
};

    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct Force
    {
        /* Typedef */
        using type = SpaceH::ProtoType<Dtype, ArraySize>;
        
        using Vector = typename type::Vector;
        
        using VectorArray = typename type::VectorArray;
        /* Typedef */
        
        void addTotal(VectorArray& acc)
        {
            const size_t size = acc.size();
            
            for(size_t i = 0 ; i < size; ++i)
            {
                acc[i] += this_acc_[i];
            }
        }
        
        const VectorArray& acc() { return this_acc_; }
        
        const Vector& acc(size_t i) { return this_acc_[i]; }
        
    protected:
        VectorArray this_acc_;
        PairForce   force_;
    };
 
    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct VelIndepForce : public Force<PairForce, Dtype, ArraySize>
    {
        /* Typedef */
        using Base = Force<PairForce, Dtype, ArraySize>;
        
        using type = typename Base::type;
        
        using Scalar = typename type::Scalar;
        
        using Vector = typename type::Vector;
        
        using ScalarArray = typename type::ScalarArray;
        
        using VectorArray = typename type::VectorArray;
        
        using IndexArray = typename type::IndexArray;
        /* Typedef */
        
        using Base::force_;
        using Base::this_acc_;
        
        constexpr static bool isVelDep{false};
        
        void calcuAcc(const ScalarArray& mass, const VectorArray& pos)
        {
            const size_t size = mass.size();
            
            for(size_t i = 0 ; i < size; ++i)
                this_acc_[i].setZero();
            
            for(size_t i = 0 ; i < size ; ++i)
            {
                for(size_t j = i + 1 ; j < size ; ++j)
                {
                    force_(this_acc_[i], this_acc_[j], mass[i], mass[j], pos[i], pos[j], pos[j] - pos[i]);
                }
            }
        }
        
        void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& chainPos, const IndexArray& chainIndex)
        {
            const size_t size = mass.size();
            
            for(size_t i = 0 ; i < size; ++i)
                this_acc_[i].setZero();
            
            for(size_t i = 0 ; i < size - 1; ++i)
            {
                force_(this_acc_[chainIndex[i]], this_acc_[chainIndex[i + 1]],
                            mass[chainIndex[i]], mass[chainIndex[i + 1]],
                            pos[chainIndex[i]], pos[chainIndex[i + 1]], chainPos[i]);

            }
            
            for(size_t i = 0 ; i < size - 2; ++i)
            {
                force_(this_acc_[chainIndex[i]], this_acc_[chainIndex[i + 2]],
                            mass[chainIndex[i]], mass[chainIndex[i + 2]],
                            pos[chainIndex[i]], pos[chainIndex[i + 2]], chainPos[i] + chainPos[i + 1]);
            }
            
            for(size_t i = 0 ; i < size ; ++i)
            {
                for(size_t j = i + 3 ; j < size; ++j)
                {
                    force_(this_acc_[chainIndex[i]], this_acc_[chainIndex[j]],
                                mass[chainIndex[i]], mass[chainIndex[j]],
                                pos[chainIndex[i]], pos[chainIndex[j]],
                                pos[chainIndex[j]] - pos[chainIndex[i]]);
                }
            }
        }
        
    };
    
    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct ExtVelIndepForce : public Force<PairForce, Dtype, ArraySize>
    {
        /* Typedef */
        using Base = Force<PairForce, Dtype, ArraySize>;
        
        using type = typename Base::type;
        
        using ScalarArray = typename type::ScalarArray;
        
        using VectorArray = typename type::VectorArray;
        
        /* Typedef */
        
        using Base::force_;
        using Base::this_acc_;
        
        constexpr static bool isVelDep{false};
        
        void calcuAcc(const ScalarArray& mass, const VectorArray& pos)
        {
            const size_t size = mass.size();
            
            for(size_t i = 0 ; i < size; ++i)
                this_acc_[i].setZero();
            
            for(size_t i = 0 ; i < size ; ++i)
            {
                force_(this_acc_[i], mass[i], pos[i]);
            }
        }
        
    };
    
    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct VelDepForce : public Force<PairForce, Dtype, ArraySize>
    {
        /* Typedef */
        using Base = Force<PairForce, Dtype, ArraySize>;
        
        using type = typename Base::type;
        
        using Scalar = typename type::Scalar;
        
        using Vector = typename type::Vector;
        
        using ScalarArray = typename type::ScalarArray;
        
        using VectorArray = typename type::VectorArray;
        
        using IndexArray = typename type::IndexArray;
        /* Typedef */
        
        using Base::force_;
        using Base::this_acc_;
        
        constexpr static bool isVelDep{true};
        
        void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
        {
            const size_t size = mass.size();
            
            for(size_t i = 0 ; i < size; ++i)
                this_acc_[i].setZero();
            
            for(size_t i = 0 ; i < size ; ++i)
            {
                for(size_t j = i + 1 ; j < size ; ++j)
                {
                    force_(this_acc_[i], this_acc_[j], mass[i], mass[j],
                                pos[i], pos[j], vel[i], vel[j],
                                pos[j] - pos[i], vel[j] - vel[i]);
                }
            }
        }
        
        void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel,
                        const VectorArray& chainPos, const VectorArray& chainVel, const IndexArray& chainIndex)
        {
            const size_t size = mass.size();
            
            for(size_t i = 0 ; i < size; ++i)
                this_acc_[i].setZero();
            
            for(size_t i = 0 ; i < size - 1; ++i)
            {
                force_(this_acc_[chainIndex[i]], this_acc_[chainIndex[i + 1]],
                            mass[chainIndex[i]], mass[chainIndex[i + 1]],
                            pos[chainIndex[i]], pos[chainIndex[i + 1]],
                            vel[chainIndex[i]], vel[chainIndex[i + 1]],
                            chainPos[i], chainVel[i]);
                
            }
            
            for(size_t i = 0 ; i < size - 2; ++i)
            {
                force_(this_acc_[chainIndex[i]], this_acc_[chainIndex[i + 2]],
                            mass[chainIndex[i]], mass[chainIndex[i + 2]],
                            pos[chainIndex[i]], pos[chainIndex[i + 2]],
                            vel[chainIndex[i]], vel[chainIndex[i + 2]],
                            chainPos[i] + chainPos[i + 1], chainVel[i] + chainVel[i + 1]);
            }
            
            for(size_t i = 0 ; i < size ; ++i)
            {
                for(size_t j = i + 3 ; j < size; ++j)
                {
                    force_(this_acc_[chainIndex[i]], this_acc_[chainIndex[j]],
                                mass[chainIndex[i]], mass[chainIndex[j]],
                                pos[chainIndex[i]], pos[chainIndex[j]],
                                vel[chainIndex[i]], vel[chainIndex[j]],
                                pos[chainIndex[j]] - pos[chainIndex[i]],
                                vel[chainIndex[j]] - vel[chainIndex[i]]);
                }
            }
        }
    };
    
    template<typename PairForce, typename Dtype, size_t ArraySize>
    struct ExtVelDepForce : public Force<PairForce, Dtype, ArraySize>
    {
        /* Typedef */
        using Base = Force<PairForce, Dtype, ArraySize>;
        
        using type = typename Base::type;
        
        using ScalarArray = typename type::ScalarArray;
        
        using VectorArray = typename type::VectorArray;
        /* Typedef */
        
        using Base::force_;
        using Base::this_acc_;
        
        constexpr static bool isVelDep{true};
        
        void calcuAcc(const ScalarArray& mass, const VectorArray& pos, const VectorArray& vel)
        {
            const size_t size = mass.size();
            
            for(size_t i = 0 ; i < size; ++i)
                this_acc_[i].setZero();
            
            for(size_t i = 0 ; i < size ; ++i)
            {
                force_(this_acc_[i], mass[i], pos[i], vel[i]);
            }
        }
        
    };
    
    template<typename Dtype, size_t ArraySize>
    struct NewtonForce
    {
        /* Typedef */
        using type = SpaceH::ProtoType<Dtype, ArraySize>;
        
        using Scalar = typename type::Scalar;
        
        using Vector = typename type::Vector;
        /* Typedef */
        inline void operator()(Vector& acc1, Vector& acc2, const Scalar m1, const Scalar m2, const Vector& pos1, const Vector& pos2, const Vector& dr)
        {
            Scalar inv_r  = dr.reNorm();
            Scalar inv_r3 = inv_r * inv_r * inv_r;
            acc1 += dr * (inv_r3 * m2);
            acc2 -= dr * (inv_r3 * m1);
        }
    };
    
}

#endif
