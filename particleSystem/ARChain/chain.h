#ifndef CHAIN_H
#define CHAIN_H
#include<vector>
#include<array>
#include<algorithm>
namespace SpaceH
{
    

namespace chain
{
/** @brief Struture to store the relative distance and index of two particles.*/
template <typename Scalar>
struct Node
{
    Scalar Rij; /**< Relative distance of two particles.*/
    size_t i; /**< Particle index.*/
    size_t j; /**< Particle index.*/
    bool   available; /**< State of node. If this node can be chained.*/
};

/** @brief Calculate the mapping index from Cartesian coordinate to chain coordinate.
 *
 *  Find the mapping index from Cartesian coordinate to chain coordinate. The chain is formed by
 *  connecting the nearest particle pairs consequently.
 *  @param pos        The array of particle position, used to calculate the distance of particle pairs.
 *  @param chainIndex The maping index needs to be calculated as a return value.
 */
template <typename VectorArray, typename IndexArray>
    void getChainIndex(const VectorArray& pos,  IndexArray& chainIndex)
    {
        using Scalar = typename VectorArray::value_type::value_type;
        
        const size_t N = pos.size();
        
        std::vector<Node<Scalar>> AdjMatrix;
        AdjMatrix.resize(N*(N - 1) / 2);
        
        createAdjMartix(pos, AdjMatrix);
        
        std::sort(AdjMatrix.begin(), AdjMatrix.end(), [&](const Node<Scalar>& Ni, const Node<Scalar>& Nj)
                  {
                      return (Ni.Rij < Nj.Rij);
                  });
        createChainIndex(AdjMatrix, chainIndex);
    }

/** @brief Create the adjoint matrix for particle pairs.
 *
 *  Create the adjoint matrix(distance of particle pairs organized by index-index matrix).
 *  @param pos        The array of particle position, used to calculate the distance of particle pairs.
 *  @param AdjMatrix  The adjoint matrix needs to be calculated as a return value.
 */
template <typename VectorArray, typename NodeArray>
void createAdjMartix(const VectorArray& pos, NodeArray& AdjMatrix )
{
    size_t N = pos.size();
    size_t k = 0;

    for(size_t i = 0 ; i < N ; ++i )
    {
        for(size_t j = i + 1 ; j < N ; ++j)
        {
            AdjMatrix[k].Rij       = (pos[i] - pos[j]).norm();
            AdjMatrix[k].i         = i;
            AdjMatrix[k].j         = j;
            AdjMatrix[k].available = true;
            k++;
        }
    }
}

/** @brief Create mapping index from adjoint matrix.
 *
 *  Create mapping index from sorted elements of adjoint matrix and connect them to a chain consequently.
 *  @param AdjMatrix  The adjoint matrix.
 *  @param chainIndex The maping index needs to be calculated as a return value.
 */
template <typename NodeArray, typename IndexArray>
void createChainIndex(NodeArray& AdjMatrix, IndexArray& chainIndex)
{
    size_t N = chainIndex.size();
    size_t chainedNumber = 0;
    size_t AdjSize = N * (N - 1) / 2;
    std::vector<size_t> Index;
    Index.reserve(N);
    Index.push_back(AdjMatrix[0].i);
    Index.push_back(AdjMatrix[0].j);
    AdjMatrix[0].available = false;
    chainedNumber++;

    for(size_t k = 0 ; k < AdjSize; ++k)
    {
        if(AdjMatrix[k].available)
        {
            size_t head = Index[0];
            size_t tail = Index[chainedNumber];

            if(AdjMatrix[k].i == head)
            {
                if( Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].j) )
                {
                    Index.insert(Index.begin(), AdjMatrix[k].j);
                    chainedNumber++;

                    if(chainedNumber < N)
                        k = 0;
                    else
                        return;
                }

                AdjMatrix[k].available = false;
                continue;
            }
            else if(AdjMatrix[k].i == tail)
            {
                if( Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].j) )
                {
                    Index.push_back(AdjMatrix[k].j);
                    chainedNumber++;

                    if(chainedNumber < N)
                        k = 0;
                    else
                        return;
                }

                AdjMatrix[k].available = false;
                continue;
            }

            if(AdjMatrix[k].j == head)
            {
                if( Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].i) )
                {
                    Index.insert(Index.begin(), AdjMatrix[k].i);
                    chainedNumber++;

                    if(chainedNumber < N)
                        k = 0;
                    else
                        return;
                }

                AdjMatrix[k].available = false;
                continue;
            }
            else if(AdjMatrix[k].j == tail)
            {
                if( Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].i) )
                {
                    Index.push_back(AdjMatrix[k].i);
                    chainedNumber++;

                    if(chainedNumber < N)
                        k = 0;
                    else
                        return;
                }

                AdjMatrix[k].available = false;
                continue;
            }
        }
    }

    for(size_t i = 0 ; i < N; ++i)
        chainIndex[i] = Index[i];
}

/** @brief Check if two mapping indexes are the same.
 *
 *  Checking the identity of two chain index mappings.
 *  @param Index1  The first index array.
 *  @param Index2  The second index array.
 *  @return boolean
 *  @note  [2,4,5,3,1] is identical to [1,3,5,4,2]
 */
template <typename IndexArray>
bool isDiff(const IndexArray& Index1, const IndexArray& Index2)
{
    const size_t N = Index1.size();
    
    if(Index1[0] == Index2[0])
    {
        for(int i = 1 ; i < N; ++i)
        {
            if(Index1[i] != Index2[i])
                return true;
        }

        return false;
    }
    else if(Index1[0] == Index2[N - 1])
    {
        for(int i = 1 ; i < N; ++i)
        {
            if(Index1[i] != Index2[N - 1 - i])
                return true;
        }

        return false;
    }
    else
        return true;
}

/** @brief Update the position chain.
 *
 *  Update the position chain. Due to the evolution, the chain index mapping could change with time,
 *  this function is used to update the position chain with old chain data.
 *  @param pos        The old chain position array needs update.
 *  @param chainIndex The old chain index mapping.
 *  @param newIndex   The new chain index mapping.
 */
template <typename VectorArray, typename IndexArray>
void updateChain(VectorArray& pos,  IndexArray& chainIndex, IndexArray& newIndex)
{
    
    size_t head    = 0;
    size_t tail    = 0;
    size_t oldhead = 0;
    size_t oldtail = 0;
    size_t size    = pos.size() - 1;
    
    typename VectorArray::value_type newPos[size];
    
    for(size_t i = 0 ; i < size; i++)
    {
        head    = newIndex[i];
        tail    = newIndex[i + 1];
        oldhead = std::find(chainIndex.begin(), chainIndex.end(), head) - chainIndex.begin();
        oldtail = std::find(chainIndex.begin(), chainIndex.end(), tail) - chainIndex.begin();
        newPos[i].setZero();

        if(oldhead < oldtail)
        {
            for(size_t j = oldhead; j < oldtail; ++j)
                newPos[i] += pos[j];
        }
        else
        {
            for(size_t j = oldtail; j < oldhead; ++j)
                newPos[i] -= pos[j];
        }
    }

    for(int i = 0 ; i < size; ++i)
        pos[i] = newPos[i];
}

/** @brief Calulate the chain data from Cartesian data and chain index mapping.
 *
 *  @param data       Data in Cartesian coordinates.
 *  @param chainData  Data need to be calculated in chain coordinates.
 *  @param chainIndex Chain index mapping.
 *  @note This function should be a inverse transformation of synCartesian().
 */
template <typename VectorArray, typename IndexArray>
void synChain(const VectorArray& data, VectorArray& chainData, const IndexArray& chainIndex)
{
    const size_t N = data.size();
    
    chainData[N - 1].setZero();
    
    for(int i = 0 ; i < N - 1; ++i)
        chainData[i] = data[chainIndex[i + 1]] - data[chainIndex[i]];
}

/** @brief Calulate the Cartesian data from chain data and chain index mapping.
 *
 *  @param chainData  Data in chain coordinates.
 *  @param data       Data need to be calculated in Cartesian coordinates.
 *  @param chainIndex Chain index mapping.
 *  @note This function should be a inverse transformation of synChain().
 */
template <typename VectorArray, typename IndexArray>
void synCartesian(const VectorArray& chainData, VectorArray& data, const IndexArray& chainIndex)
{
    const size_t N = data.size();
    
    data[chainIndex[0]].setZero();

    for(int i = 1 ; i < N ; ++i)
        data[chainIndex[i]] = data[chainIndex[i - 1]] + chainData[i - 1];
}
}
}
#endif

