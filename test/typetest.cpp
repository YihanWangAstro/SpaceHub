
#include<iostream>
#include<array>

template<typename T>
struct a
{
    template<typename U, size_t N>
    using Container = std::array<U, N>;
    typedef Container<T,2> real;
};

template<typename T>
struct b
{
    template<typename U, size_t N>
    using Container = typename T::template Container<U,N>;
    typedef typename T::real real;
    void foo()
    {
      //  Container<int,2> arry;
        real tmp;
        //std::cout << arry[0] << std::endl;
    }
};

template<typename in>
struct p
{
    typename in::real bar()
    {
        std::cout << "g" << std::endl;
        return 0;
    }
};

int main()
{
    
    
    b<a<int> > s;
    s.foo();
    return 0;
}
