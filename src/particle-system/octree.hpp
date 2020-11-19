//
// Created by 王艺涵 on 4/15/19.
//

#pragma once

#include <memory>
namespace space::octree {

    template <typename T>
    class Node {
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

    enum class region { flu, fld, fru, frd, blu, bld, bru, brd };

    template <typename T>
    class Octree {
       public:
        using value_type = T;

        template <typename STL>
        Octree(STL& stl) {
            for (auto& d : stl) {
                insert();
            }
        }

       private:
        void insert(std::unique_ptr<Node<T>>& node, T const& data) {
            if (node == nullptr) {
                node = std::make_unique<Node<T>>();
                node.data = &data;
            } else {
                auto loc = get_region(node.data->pos, data.pos);
                switch (loc) {
                    case region::flu:
                        insert(node.flu, data);
                        break;
                    case region::fld:
                        insert(node.fld, data);
                        break;
                    case region::fru:
                        insert(node.fru, data);
                        break;
                    case region::frd:
                        insert(node.frd, data);
                        break;
                    case region::blu:
                        insert(node.blu, data);
                        break;
                    case region::bld:
                        insert(node.bld, data);
                        break;
                    case region::bru:
                        insert(node.bru, data);
                        break;
                    case region::brd:
                        insert(node.brd, data);
                        break;
                }
            }
        }

        template <typename Vector>
        region get_region(Vector const& loc, Vector const& p) {
            if (p.x > loc.x) {
                if (p.y > loc.y) {
                    if (p.z > loc.z) {
                        return region::fru;
                    } else {
                        return region::frd;
                    }
                } else {
                    if (p.z > loc.z) {
                        return region::flu;
                    } else {
                        return region::fld;
                    }
                }
            } else {
                if (p.y > loc.y) {
                    if (p.z > loc.z) {
                        return region::bru;
                    } else {
                        return region::brd;
                    }
                } else {
                    if (p.z > loc.z) {
                        return region::blu;
                    } else {
                        return region::bld;
                    }
                }
            }
        }

        std::unique_ptr<Node<T>> root{nullptr};
    };
}  // namespace space::octree
