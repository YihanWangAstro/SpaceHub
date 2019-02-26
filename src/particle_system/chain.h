#ifndef CHAIN_H
#define CHAIN_H

#include<vector>
#include<array>
#include<algorithm>

namespace SpaceH {


    namespace chain {
/** @brief Struture to store the relative distance and index of two particles.*/
        template<typename Scalar>
        struct Node {
            Node(Scalar d, size_t ix, size_t jx, bool a) :
            Rij(d), i(ix), j(jx), available(a){}
            Node() = default;
            Scalar Rij; /**< Relative distance of two particles.*/
            size_t i; /**< Particle index.*/
            size_t j; /**< Particle index.*/
            bool available; /**< State of node. If this node can be chained.*/
        };

/** @brief Create the adjoint matrix for particle pairs.
 *
 *  Create the adjoint matrix(distance of particle pairs organized by index-index matrix).
 *  @param pos        The array of particle position, used to calculate the distance of particle pairs.
 *  @param AdjMatrix  The adjoint matrix needs to be calculated as a return value.
 */
        template<typename VectorArray, typename NodeArray>
        void createAdjMartix(const VectorArray &pos, NodeArray &AdjMatrix) {
            using Scalar = typename VectorArray::value_type::value_type;
            size_t N = pos.size();
            size_t k = 0;
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = i + 1; j < N; ++j) {
                    AdjMatrix[k++] = Node<Scalar>(distance(pos[i],pos[j]), i, j, true);
                }
            }
        }

        /** @brief Create mapping index from adjoint matrix.
     *
     *  Create mapping index from sorted elements of adjoint matrix and connect them to a chain consequently.
     *  @param AdjMatrix  The adjoint matrix.
     *  @param chainIndex The maping index needs to be calculated as a return value.
     */
        template<typename NodeArray, typename IndexArray>
        void createChainIndex(NodeArray &AdjMatrix, IndexArray &chainIndex) {
            size_t N = chainIndex.size();
            size_t chainedNumber = 0;
            size_t AdjSize = N * (N - 1) / 2;
            std::vector<size_t> Index;
            Index.reserve(N);
            Index.emplace_back(AdjMatrix[0].i);
            Index.emplace_back(AdjMatrix[0].j);
            AdjMatrix[0].available = false;
            chainedNumber++;

            for (size_t k = 0; k < AdjSize; ++k) {
                if (AdjMatrix[k].available) {
                    size_t head = Index[0];
                    size_t tail = Index[chainedNumber];

                    if (AdjMatrix[k].i == head) {
                        if (Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].j)) {
                            Index.insert(Index.begin(), AdjMatrix[k].j);
                            chainedNumber++;

                            if (chainedNumber < N)
                                k = 0;
                            else
                                return;
                        }

                        AdjMatrix[k].available = false;
                        continue;
                    } else if (AdjMatrix[k].i == tail) {
                        if (Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].j)) {
                            Index.emplace_back(AdjMatrix[k].j);
                            chainedNumber++;

                            if (chainedNumber < N)
                                k = 0;
                            else
                                return;
                        }

                        AdjMatrix[k].available = false;
                        continue;
                    }

                    if (AdjMatrix[k].j == head) {
                        if (Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].i)) {
                            Index.insert(Index.begin(), AdjMatrix[k].i);
                            chainedNumber++;

                            if (chainedNumber < N)
                                k = 0;
                            else
                                return;
                        }

                        AdjMatrix[k].available = false;
                        continue;
                    } else if (AdjMatrix[k].j == tail) {
                        if (Index.end() == std::find(Index.begin(), Index.end(), AdjMatrix[k].i)) {
                            Index.emplace_back(AdjMatrix[k].i);
                            chainedNumber++;

                            if (chainedNumber < N)
                                k = 0;
                            else
                                return;
                        }

                        AdjMatrix[k].available = false;
                        continue;
                    }
                }
            }

            for (size_t i = 0; i < N; ++i)
                chainIndex[i] = Index[i];
        }

/** @brief Calculate the mapping index from Cartesian coordinate to chain coordinate.
 *
 *  Find the mapping index from Cartesian coordinate to chain coordinate. The chain is formed by
 *  connecting the nearest particle pairs consequently.
 *  @param pos        The array of particle position, used to calculate the distance of particle pairs.
 *  @param chainIndex The maping index needs to be calculated as a return value.
 */
        template<typename VectorArray, typename IndexArray>
        void getChainIndex(const VectorArray &pos, IndexArray &chainIndex) {
            using Scalar = typename VectorArray::value_type::value_type;

            const size_t N = pos.size();
            if(chainIndex.size() != N) {
                chainIndex.resize(pos.size());
            }
            std::vector<Node<Scalar>> AdjMatrix;
            AdjMatrix.resize(N * (N - 1) / 2);

            createAdjMartix(pos, AdjMatrix);

            std::sort(AdjMatrix.begin(), AdjMatrix.end(), [&](const Node<Scalar> &Ni, const Node<Scalar> &Nj) {
                return (Ni.Rij < Nj.Rij);
            });
            createChainIndex(AdjMatrix, chainIndex);
        }


/** @brief Check if two mapping indexes are the same.
 *
 *  Checking the identity of two chain index mappings.
 *  @param index1  The first index array.
 *  @param index2  The second index array.
 *  @return boolean
 *  @note  [2,4,5,3,1] is identical to [1,3,5,4,2]
 */
        template<typename IndexArray>
        bool is_diff(const IndexArray &index1, const IndexArray &index2) {
            const size_t size = index1.size();

            if (index1[0] == index2[0]) {
                for (int i = 1; i < size; ++i) {
                    if (index1[i] != index2[i])
                        return true;
                }

                return false;
            } else if (index1[0] == index2[size - 1]) {
                for (int i = 1; i < size; ++i) {
                    if (index1[i] != index2[size - 1 - i])
                        return true;
                }

                return false;
            } else
                return true;
        }

/** @brief Update the position chain.
 *
 *  Update the position chain. Due to the evolution, the chain index mapping could change with time,
 *  this function is used to update the position chain with old chain data.
 *  @param chain        The old chain position array needs update.
 *  @param index The old chain index mapping.
 *  @param newIndex   The new chain index mapping.
 */
        template<typename DataArray, typename IndexArray>
        void update_chain(DataArray &chain, IndexArray &index, IndexArray &newIndex) {
            size_t size = chain.size() - 1;
            typename DataArray::value_type newPos[size];

            size_t head0 = std::find(index.begin(), index.end(), newIndex[0]) - index.begin();
            typename DataArray::value_type headPos = chain[size];
            for (int i = 0; i < head0; ++i)
                headPos += chain[i];
            chain[size] = headPos;
            //pos[size].setZero();

            for (size_t i = 0; i < size; i++) {
                size_t head = newIndex[i];
                size_t tail = newIndex[i + 1];
                size_t old_head = std::find(index.begin(), index.end(), head) - index.begin();
                size_t old_tail = std::find(index.begin(), index.end(), tail) - index.begin();
                newPos[i] = 0;

                if (old_head < old_tail) {
                    for (size_t j = old_head; j < old_tail; ++j)
                        newPos[i] += chain[j];
                } else {
                    for (size_t j = old_tail; j < old_head; ++j)
                        newPos[i] -= chain[j];
                }
            }

            for (int i = 0; i < size; ++i)
                chain[i] = newPos[i];
        }


        template<typename DataArray, typename IndexArray>
        void to_chain_coord(DataArray const &data, DataArray &chain, IndexArray const &index) {
            const size_t size = data.size();
            //chain[N - 1] = 0;
            chain[size-1] = data[index[0]];
            for (size_t i = 0; i < size - 1; ++i)
                chain[i] = data[index[i + 1]] - data[index[i]];
        }

        template<typename DataArray, typename IndexArray>
        void to_cartesian_coord(DataArray const &chain, DataArray &data, IndexArray const &index) {
            const size_t size = data.size();

            //data[index[0]] = 0;
            data[index[0]] = chain[size-1];
            for (size_t i = 1; i < size; ++i)
                data[index[i]] = data[index[i - 1]] + chain[i - 1];
        }
    }
}
#endif

