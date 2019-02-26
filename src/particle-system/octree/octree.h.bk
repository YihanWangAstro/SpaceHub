#ifndef OCTREE_H
#define OCTREE_H
namespace SpaceH
{
    
    namespace OcTree
    {
        template<typename T>
        struct Oct
        {
            using value_type = T;
            
            T trf{static_cast<T>(0)};
            T trb{static_cast<T>(0)};
            T tlf{static_cast<T>(0)};
            T tlb{static_cast<T>(0)};
            T brf{static_cast<T>(0)};
            T brb{static_cast<T>(0)};
            T blf{static_cast<T>(0)};
            T blb{static_cast<T>(0)};
        };
        
        template<typename T>
        struct Coord
        {
            using value_type = T;
            
            T xmin{0};
            T xmax{0};
            T ymin{0};
            T ymax{0};
            T zmin{0};
            T zmax{0};
        };
        
        template<typename T>
        struct ListNode
        {
            T data;
            ListNode* next{nullptr};
        };
        
        template<typename T>
        struct List
        {
            using Node = ListNode<T>;
            
            size_t size{0};
            Node *head{nullptr};
            Node *tail{nullptr};
        };
        template<typename Dtype, typename CoordType>
        struct OctreeNode
        {
            Node(size_t num, Coord<CoordType>& coords) : number(num), coord(coords){}
            
            Dtype*   data{nullptr};
            Coord<CoordType> coord;
            Oct<Node*>    children;
            size_t       number{0};
            bool     isLeaf{false};
        };
        
        template<typename T>
        struct emptyParentNode
        {
            T* operator()(T* dataSet){ return nullptr; }
        };
        
        template<typename List, typename T>
        void splitOct(List& rawList, Oct<List>& Lists, T xmid, T ymid, T zmid)
        {
            Oct<ListNode*> tail;
            ListNode* iter = rawList;
            for(;iter != nullptr;)
            {
                if(iter->data.z() > zmid && iter->data.y() > ymid && iter->data.x() > xmid)
                {
                    sublist.trf.add(iter);
                }else if(iter->data.z() > zmid && iter->data.y() > ymid && iter->data.x() <= xmid)
                {
                    sublist.trb.add(iter);
                }else if(iter->data.z() > zmid && iter->data.y() <= ymid && iter->data.x() > xmid)
                {
                    sublist.tlf.add(iter);
                }else if(iter->data.z() > zmid && iter->data.y() <= ymid && iter->data.x() <= xmid)
                {
                    sublist.tlb.add(iter);
                }else if(iter->data.z() <= zmid && iter->data.y() > ymid && iter->data.x() > xmid)
                {
                    sublist.brf.add(iter);
                }else if(iter->data.z() <= zmid && iter->data.y() > ymid && iter->data.x() <= xmid)
                {
                    sublist.brb.add(iter);
                }else if(iter->data.z() <= zmid && iter->data.y() <= ymid && iter->data.x() > xmid)
                {
                    sublist.blf.add(iter);
                }else if(iter->data.z() <= zmid && iter->data.y() <= ymid && iter->data.x() <= xmid)
                {
                    sublist.blb.add(iter);
                }else
                {
                    SpaceH::errMsg("unexpected coordinates in octree!", __FILE__, __LINE__);
                }
                iter = iter->next;
            }
                
        };
        
        template<typename Dtype, typename CoordType, typename ParentFun = emptyParentNode<Dtype>>
        struct Octree
        {
            using TreeNode = OctreeNode<Dtype, CoordType>;
            using Box      = Coord<CorrdType>;
            using List = List<Dtype>;
        private:
            void create(TreeNode* &parent, List &list, CoordType xmin, CoordType xmax, CoordType ymin, CoordType ymax, CoordType zmin, CoordType zmax)
            {
                if(list.number != 0)
                {
                    parent = new TreeNode;
                    parent->number = list.size();
                    parent->coord.xmin = xmin;
                    parent->coord.xmax = xmax;
                    parent->coord.ymin = ymin;
                    parent->coord.ymax = ymax;
                    parent->coord.zmin = zmin;
                    parent->coord.zmax = zmax;
                    
                    if(parent->number == 1)
                    {
                        parent->data = list.front();
                        parent->isleaf = true;
                    }
                    else
                    {
                        parent->data = getParent(parent);
                        
                        CoordType xmid = 0.5*(xmin + xmax);
                        CoordType ymid = 0.5*(ymin + ymax);
                        CoordType zmid = 0.5*(zmin + zmax);
                        
                        Oct<List> sublist;
                        
                        splitOct(list, sublist, xmid, ymid, zmid);
                        
                        create(parent->children.trf, sublist.trf, xmid, xmax, ymid, ymax, zmid, zmax);
                        create(parent->children.trb, sublist.trb, xmin, xmid, ymid, ymax, zmid, zmax);
                        create(parent->children.tlf, sublist.tlf, xmid, xmax, ymin, ymid, zmid, zmax);
                        create(parent->children.tlb, sublist.tlb, xmin, xmid, ymin, ymid, zmid, zmax);
                        create(parent->children.brf, sublist.brf, xmid, xmax, ymid, ymax, zmin, zmid);
                        create(parent->children.brb, sublist.brb, xmin, xmid, ymid, ymax, zmin, zmid);
                        create(parent->children.blf, sublist.blf, xmid, xmax, ymin, ymid, zmin, zmid);
                        create(parent->children.blb, sublist.blb, xmin, xmid, ymin, ymid, zmin, zmid);
                    }
                }
                else
                {
                    parent = nullptr;
                }
            }
            
            void classify(Oct<List>& sublist, ListNode* & list, CoordType xmid, CoordType ymid, CoordType zmid)
            {
                typename List::NodeType *iter = list->next;
                
            }
        private:
            TreeNode* root{nullptr};
            ParentFun getParent;
        }
    }
}
#endif
}
