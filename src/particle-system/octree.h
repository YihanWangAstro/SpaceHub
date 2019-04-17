//
// Created by 王艺涵 on 4/15/19.
//

#ifndef SPACEHUB_OCTREE_H
#define SPACEHUB_OCTREE_H

#include <memory>
namespace space::octree{

    template<typename T>
    class Node{
        T* data{nullptr};
        std::unique_ptr<Node> flu{nullptr};
        std::unique_ptr<Node> fld{nullptr};
        std::unique_ptr<Node> fru{nullptr};
        std::unique_ptr<Node> frd{nullptr};
        std::unique_ptr<Node> blu{nullptr};
        std::unique_ptr<Node> bld{nullptr};
        std::unique_ptr<Node> bru{nullptr};
        std::unique_ptr<Node> brd{nullptr};
    };

    enum class region{
        flu, fld, fru, frd, blu, bld, bru, brd
    };

    template<typename T>
    class Octree{
    public:
        using value_type = T;

        template <typename STL>
        Octree(STL& stl) {
            for(auto& d : stl) {
                insert()
            }
        }
    private:
        void insert(std::unique_ptr<Node<T>>& node, T const& data){
            if(node == nullptr) {
                node = new Node<T>;
                node.data = &data;
            } else {
                auto loc = get_region(*(node.data), data);
                switch (loc){
                    case region::flu :
                        insert(node.flu, data);
                        break;
                    case region::fld :
                        insert(node.fld, data);
                        break;
                    case region::fru :
                        insert(node.fru, data);
                        break;
                    case region::frd :
                        insert(node.frd, data);
                        break;
                    case region::blu :
                        insert(node.blu, data);
                        break;
                    case region::bld :
                        insert(node.bld, data);
                        break;
                    case region::bru :
                        insert(node.bru, data);
                        break;
                    case region::brd :
                        insert(node.brd, data);
                        break;
                }
            }
        }
        region get_region(T const& loc, T const& data){

        }
        std::unique_ptr<Node<T>> root{nullptr};
    };
}
#endif //SPACEHUB_OCTREE_H
